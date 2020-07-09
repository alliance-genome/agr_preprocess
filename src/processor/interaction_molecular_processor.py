import logging
import multiprocessing
import csv
import re
import json
import os
import sys
import urllib.request
import requests
from tqdm import tqdm
from datetime import datetime
from string import Template


from processor import Processor

logger = logging.getLogger(__name__)


class HeaderTemplate(Template):
    delimiter = '%'


class InteractionMolecularProcessor(Processor):

    def __init__(self, configs):
        super().__init__()
        self.data_type_configs = configs
        self.master_gene_set = set()
        self.master_crossreference_dictionary = dict()
        self.master_crossreference_dictionary['UniProtKB'] = dict()
        self.master_crossreference_dictionary['ENSEMBL'] = dict()
        self.master_crossreference_dictionary['NCBI_Gene'] = dict()
        self.master_crossreference_dictionary['RefSeq'] = dict()
        self.biogrid_rna_set = set()
        self.biogrid_doi_dict = dict()
        self.coronavirus_injection_dict = dict()
        self.output_dir = '/usr/src/app/output/'
        self.download_dir = '/usr/src/app/download_molecular/'


    def _load_and_process_data(self):
        logger.info("in InteractionMolecularProcessor")

        source_filepaths = dict()
        interaction_source_config = self.data_type_configs[0]
        for sub_type in interaction_source_config.get_sub_type_objects():
            sub_type_name = sub_type.get_sub_data_type()
            sub_type_filepath = sub_type.get_filepath()
            source_filepaths[sub_type_name] = sub_type_filepath

        for sub_type in source_filepaths:
            logger.debug("Source subtype %s filepath %s" % (sub_type, source_filepaths[sub_type]))

        bgi_filepaths = dict()
        interaction_source_config = self.data_type_configs[1]
        for sub_type in interaction_source_config.get_sub_type_objects():
            sub_type_name = sub_type.get_sub_data_type()
            sub_type_filepath = sub_type.get_filepath()
            bgi_filepaths[sub_type_name] = sub_type_filepath

        for sub_type in bgi_filepaths:
            logger.debug("BGI subtype %s filepath %s" % (sub_type, bgi_filepaths[sub_type]))

        interactions_molecular = InteractionMolecularProcessor(self.data_type_configs)
        interactions_molecular.parse_bgi_json()
        interactions_molecular.get_data()
        interactions_molecular.validate_and_upload_files_to_fms()


    def parse_bgi_json(self):
        # We're populating a rather large dictionary to use for looking up Alliance genes by their crossreferences.
        # Edit the list below if you'd like to add more crossreferences to the dictionary.
        # The key of the dictionary is the crossreference and the value is the Alliance gene to which it resolves.
        #
        # We're also populating the "master gene set" for gene lookups later.
        logger.info('Populating master gene set and crossreferences from JSON.')

        bgi_filepaths = dict()
        interaction_source_config = self.data_type_configs[1]
        for sub_type in interaction_source_config.get_sub_type_objects():
            sub_type_name = sub_type.get_sub_data_type()
            sub_type_filepath = sub_type.get_filepath()
            bgi_filepaths[sub_type_name] = sub_type_filepath

        for sub_type in bgi_filepaths:
            logger.info("BGI subtype %s filepath %s" % (sub_type, bgi_filepaths[sub_type]))
            filepath = bgi_filepaths[sub_type]
            with open(filepath) as json_file:
                data = json.load(json_file)
                logger.info('Scanning {}'.format(filepath))
                for item in tqdm(data['data']):
                    gene_identifier = item['basicGeneticEntity']['primaryId']
                    self.master_gene_set.add(gene_identifier)

                    for xref in item['basicGeneticEntity']['crossReferences']:
                        cross_ref_record = None
                        cross_ref_prefix = None
                        if xref['id'].startswith('NCBI_Gene'):
                            # Modify the cross reference ID to match the PSI MITAB format if necessary.
                            # So far, this is just converting 'NCBI_Gene' to 'entrez gene/locuslink'.
                            cross_ref_prefix = 'NCBI_Gene'
                            cross_ref_record_split = xref['id'].split(':')[1]
                            cross_ref_record = 'entrez gene/locuslink:' + cross_ref_record_split
                        elif xref['id'].startswith('UniProtKB'):
                            cross_ref_prefix = 'UniProtKB'
                            cross_ref_record = xref['id']
                        elif xref['id'].startswith('ENSEMBL'):
                            cross_ref_prefix = 'ENSEMBL'
                            cross_ref_record = xref['id']
                        elif xref['id'].startswith('RefSeq'):
                            cross_ref_prefix = 'RefSeq'
                            cross_ref_record = xref['id']

                        # The crossreference dictionary is a list of genes linked to a single crossreference.
                        # Append the gene if the crossref dict entry exists.
                        # Otherwise, create a list and append the entry.
                        if cross_ref_record is not None:
                            if cross_ref_record.lower() in self.master_crossreference_dictionary[cross_ref_prefix]:
                                self.master_crossreference_dictionary[cross_ref_prefix][cross_ref_record.lower()]\
                                    .append(gene_identifier)
                            else:
                                self.master_crossreference_dictionary[cross_ref_prefix][cross_ref_record.lower()] = []
                                self.master_crossreference_dictionary[cross_ref_prefix][cross_ref_record.lower()].append(
                                    gene_identifier)

        # hardcoded exception that does not have a resolvable crossReference identifier in BGI
        self.master_crossreference_dictionary['NCBI_Gene']['entrez gene/locuslink:ORF9c'.lower()] = ['RefSeq:P0DTD3']

                             # The ids in PSI-MITAB files are lower case, hence the .lower() used above.
        logger.info('Done.')

    def resolve_identifiers_by_row(self, row, mapped_out):
        interactor_A_rows = [0, 2, 4, 22]
        interactor_B_rows = [1, 3, 5, 23]

        interactor_A_resolved = False
        interactor_B_resolved = False

        taxon1 = re.search(r'taxid:\d+', row[9]).group(0)
        taxon2 = re.search(r'taxid:\d+', row[10]).group(0)

        interactor_A_list = []
        for row_entry in interactor_A_rows:
            if '|' in row[row_entry]:
                interactors_A = row[row_entry].split('|')
                for interactor_A in interactors_A:
                    interactor_A = interactor_A.replace("uniprot/swiss-prot:", "uniprotkb:")
                    interactor_A_list.append(interactor_A)
            else:
                interactor_A_list.append(row[row_entry])

        interactor_B_list = []
        for row_entry in interactor_B_rows:
            if '|' in row[row_entry]:
                interactors_B = row[row_entry].split('|')
                for interactor_B in interactors_B:
                    interactor_B = interactor_B.replace("uniprot/swiss-prot:", "uniprotkb:")
                    interactor_B_list.append(interactor_B)
            else:
                interactor_B_list.append(row[row_entry])

        interactor_A_resolved, A_resolved_id = self.resolve_identifier(interactor_A_list, taxon1)
        interactor_B_resolved, B_resolved_id = self.resolve_identifier(interactor_B_list, taxon2)

        if A_resolved_id is not None and B_resolved_id is not None:
            mapped_output_rows = [row[13], A_resolved_id, B_resolved_id]
            mapped_out.writerow(mapped_output_rows)
            if taxon1 == 'taxid:2697049':	# for sars-cov-2 first match on refseq
                if A_resolved_id in self.coronavirus_injection_dict:
                    row[0] = self.coronavirus_injection_dict[A_resolved_id][0]
                    row[2] = self.coronavirus_injection_dict[A_resolved_id][1]
                    row[4] = self.coronavirus_injection_dict[A_resolved_id][2]
            if taxon2 == 'taxid:2697049':	# for sars-cov-2 first match on refseq
                if B_resolved_id in self.coronavirus_injection_dict:
                    row[1] = self.coronavirus_injection_dict[B_resolved_id][0]
                    row[3] = self.coronavirus_injection_dict[B_resolved_id][1]
                    row[5] = self.coronavirus_injection_dict[B_resolved_id][2]

        return interactor_A_resolved, interactor_B_resolved, row


    def resolve_identifier(self, interactor_list, taxon):

        # The order of this list is important.
        list_of_crossref_regex_to_search = [
            'uniprotkb:[\\w\\d_-]*$',
            'ensembl:[\\w\\d_-]*$',
            'entrez gene/locuslink:.*$'
        ]
        if taxon == 'taxid:2697049':	# for sars-cov-2, first match on refseq
            list_of_crossref_regex_to_search.insert(0, 'refseq:[\\w\\d_-]*$')

        # Prioritize matching based on regex priority instead of column order
        for regex_entry in list_of_crossref_regex_to_search:
            for individual_entry in interactor_list:
                # For use in wormbase / flybase lookups.
                # If we run into an IndexError, there's no identifier to resolve and we return False.
                try:
                    entry_stripped = individual_entry.split(':')[1]
                except IndexError:
                    continue

                # uniprotkb: could have trailing '-<something>' that should be stripped
                if individual_entry.startswith('uniprotkb:'):
                    if taxon != 'taxid:2697049':		# if not sars-cov-2 strip hyphen and after
                        individual_entry = individual_entry.split('-')[0]

                prefixed_identifier = None
                logger.debug('resolving individual_entry : %s ; taxon : %s' % (individual_entry, taxon))

                if entry_stripped.startswith('WB'):
                    prefixed_identifier = 'WB:' + entry_stripped
                    if prefixed_identifier in self.master_gene_set:
                        return True, prefixed_identifier
                    else:
                        logger.debug('resolved WB False : ' + prefixed_identifier)
                        return False, None
                elif entry_stripped.startswith('FB'):
                    prefixed_identifier = 'FB:' + entry_stripped
                    if prefixed_identifier in self.master_gene_set:
                        logger.debug('resolved FB False : ' + prefixed_identifier)
                        return True, prefixed_identifier
                    else:
                        return False, None

                regex_output = re.findall(regex_entry, individual_entry)
                if regex_output is not None:
                    for regex_match in regex_output: # We might have multiple regex matches. Search them all against our crossreferences.
                        identifier = regex_match
                        for crossreference_type in self.master_crossreference_dictionary.keys():
                            # Using lowercase in the identifier to be consistent with Alliance lowercase identifiers.
                            if identifier.lower() in self.master_crossreference_dictionary[crossreference_type]:
                                if taxon == 'taxid:2697049':	# for sars-cov-2 get the first primaryId value that matches the identifier in the BGI
                                    return True, self.master_crossreference_dictionary[crossreference_type][identifier.lower()][0].lower()	# Return 'True' if we find an entry.
                                else:				# for other taxons return the identifier that matches in the BGI
                                    return True, identifier.lower() # Return 'True' if we find an entry.

        # If we can't resolve any of the crossReferences, return None
        interactor_list_string = "\t".join(interactor_list)
        logger.debug('resolved default False : ' + interactor_list_string)
        return False, None


    def unzip_to_filename(self, filename_zip, filename, datasubtype):
        logger.info('Extracting file {} with unzip into {}'.format(filename_zip, filename))
        os.system('unzip -o {} -d {}tmp/'.format(filename_zip, self.download_dir))
        logger.info('Renaming extracted file.')
        if datasubtype == 'IMEX':  # Special exception for IMEX because it's 2 files.
            os.system('mv {}tmp/intact.txt {}'.format(self.download_dir, filename))
            os.system('rm {}tmp/*'.format(self.download_dir))
        else:
            os.system('mv {}tmp/* {}'.format(self.download_dir, filename))

    def read_coronavirus_injection(self, coronavirus_injection_filename):
        with open(coronavirus_injection_filename, 'r', encoding='utf-8') as coronavirus_injection_in:
            csv_reader = csv.reader(coronavirus_injection_in, delimiter='\t', quoting=csv.QUOTE_NONE)
            next(csv_reader, None) # Skip the headers
            for row in tqdm(csv_reader):
                if row[0].lower().startswith("refseq"):
                    logger.debug('Mapping coronavirus_injection_dict %s' % (row[0].lower()))
                    self.coronavirus_injection_dict[row[0].lower()] = [ row[0], row[1], row[2] ]

    def create_biogrid_mappings(self, tab30_filename):
        with open(tab30_filename, 'r', encoding='utf-8') as tab30in:
            csv_reader = csv.reader(tab30in, delimiter='\t', quoting=csv.QUOTE_NONE)
            next(csv_reader, None) # Skip the headers
            for row in tqdm(csv_reader):
                if 'RNA' in row[11]:
                    self.biogrid_rna_set.add(row[0])
                if row[14].lower().startswith("doi"):
                    biogrid_key = 'biogrid:' + row[0]
                    self.biogrid_doi_dict[biogrid_key] = row[14]


    def get_data(self):
        source_filepaths = dict()
        interaction_source_config = self.data_type_configs[0]
        for sub_type in interaction_source_config.get_sub_type_objects():
            sub_type_name = sub_type.get_sub_data_type()
            sub_type_filepath = sub_type.get_filepath()
            source_filepaths[sub_type_name] = sub_type_filepath

        for sub_type in source_filepaths:
            logger.info("Source subtype %s filepath %s" % (sub_type, source_filepaths[sub_type]))

        wormbase_filename = source_filepaths['WB-MOL']
        flybase_filename = source_filepaths['FB-MOL']
        imex_filename_zip = source_filepaths['IMEX']
        imex_filename = self.download_dir + 'INTERACTION-MOL_IMEX'
        biogrid_filename_zip = source_filepaths['BIOGRID']
        biogrid_filename = self.download_dir + 'INTERACTION-MOL_BIOGRID'
        tab30_filename_zip = source_filepaths['BIOGRID-TAB']
        tab30_filename = self.download_dir + 'INTERACTION-MOL_BIOGRID-TAB'

        # comment this out on runs after the first one, to save time on unzipping large files
        self.unzip_to_filename(imex_filename_zip, imex_filename, 'IMEX')
        self.unzip_to_filename(biogrid_filename_zip, biogrid_filename, 'BIOGRID')
        self.unzip_to_filename(tab30_filename_zip, tab30_filename, 'BIOGRID-TAB')

        # To use smaller subset files
        # imex_filename = self.output_dir + 'IMEx_SARS-CoV-2-human-only_interactions.mitab.txt'
        # biogrid_filename = self.output_dir + 'BioGRID_SARS-CoV-2-human-only_interactions.mitab.txt'

        self.create_biogrid_mappings(tab30_filename)

        coronavirus_injection_filename = 'curation/coronavirus_injection.tsv'
        self.read_coronavirus_injection(coronavirus_injection_filename)

        # The order of this list is important.
        parsing_list = [wormbase_filename, flybase_filename, biogrid_filename, imex_filename]

        taxon_species_set = (
            'taxid:10116',
            'taxid:9606',
            'taxid:10090',
            'taxid:6239',
            'taxid:559292',
            'taxid:7955',
            'taxid:7227',
            'taxid:4932',
            'taxid:307796',
            'taxid:643680',
            'taxid:574961',
            'taxid:285006',
            'taxid:545124',
            'taxid:2697049',
            'taxid:764097',
            '-')
        possible_yeast_taxon_set = ('taxid:4932', 'taxid:307796', 'taxid:643680', 'taxid:574961', 'taxid:285006', 'taxid:545124', 'taxid:764097')
        interaction_exclusion_set = ('psi-mi:\"MI:0208\"', 'psi-mi:\"MI:0794\"', 'psi-mi:\"MI:0796\"', 'psi-mi:\"MI:0799\"')
        interactor_type_exclusion_set = ('psi-mi:\"MI:0328\"', 'psi-mi:\"MI:1302\"', 'psi-mi:\"MI:1304\"', 'psi-mi:\"MI:0680\"')

        # Load genes and uniprot csv files retrieved from the Alliance.

        psi_mi_tab_header = [
            '#ID(s) interactor A',
            'ID(s) interactor B',
            'Alt. ID(s) interactor A',
            'Alt. ID(s) interactor B',
            'Alias(es) interactor A',
            'Alias(es) interactor B',
            'Interaction detection method(s)',
            'Publication 1st author(s)',
            'Publication Identifier(s)',
            'Taxid interactor A',
            'Taxid interactor B',
            'Interaction type(s)',
            'Source database(s)',
            'Interaction identifier(s)',
            'Confidence value(s)',
            'Expansion method(s)',
            'Biological role(s) interactor A',
            'Biological role(s) interactor B',
            'Experimental role(s) interactor A',
            'Experimental role(s) interactor B',
            'Type(s) interactor A',
            'Type(s) interactor B',
            'Xref(s) interactor A',
            'Xref(s) interactor B',
            'Interaction Xref(s)',
            'Annotation(s) interactor A',
            'Annotation(s) interactor B',
            'Interaction annotation(s)',
            'Host organism(s)',
            'Interaction parameter(s)',
            'Creation date',
            'Update date',
            'Checksum(s) interactor A',
            'Checksum(s) interactor B',
            'Interaction Checksum(s)',
            'Negative',
            'Feature(s) interactor A',
            'Feature(s) interactor B',
            'Stoichiometry(s) interactor A',
            'Stoichiometry(s) interactor B',
            'Identification method participant A',
            'Identification method participant B'
        ]

        publication_tracking_dict = {}

        # Open all of the output files.
        with open(self.output_dir + 'alliance_molecular_interactions.tsv', 'w', encoding='utf-8') as tsvout, \
             open(self.output_dir + 'alliance_molecular_interactions_sarscov2.tsv', 'w', encoding='utf-8') as sarscov2_out, \
             open(self.output_dir + 'alliance_molecular_interactions_fly.tsv', 'w', encoding='utf-8') as fb_out, \
             open(self.output_dir + 'alliance_molecular_interactions_worm.tsv', 'w', encoding='utf-8') as wb_out, \
             open(self.output_dir + 'alliance_molecular_interactions_zebrafish.tsv', 'w', encoding='utf-8') as zfin_out, \
             open(self.output_dir + 'alliance_molecular_interactions_yeast.tsv', 'w', encoding='utf-8') as sgd_out, \
             open(self.output_dir + 'alliance_molecular_interactions_rat.tsv', 'w', encoding='utf-8') as rgd_out, \
             open(self.output_dir + 'alliance_molecular_interactions_mouse.tsv', 'w', encoding='utf-8') as mgi_out, \
             open(self.output_dir + 'alliance_molecular_interactions_human.tsv', 'w', encoding='utf-8') as human_out, \
             open(self.output_dir + 'molecular_interactions_skipped_entries.txt', 'w', encoding='utf-8') as skipped_out, \
             open(self.output_dir + 'molecular_interactions_mapped_entries.txt', 'a+', encoding='utf-8') as mapped_out:

            mapped_out = csv.writer(mapped_out, quotechar = '', quoting=csv.QUOTE_NONE, delimiter='\t')
            tsvout = csv.writer(tsvout, quotechar = '', quoting=csv.QUOTE_NONE, delimiter='\t')
            skipped_out = csv.writer(skipped_out, quotechar = '', quoting=csv.QUOTE_NONE, delimiter='\t')
            sarscov2_out = csv.writer(sarscov2_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            fb_out = csv.writer(fb_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            wb_out = csv.writer(wb_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            zfin_out = csv.writer(zfin_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            sgd_out = csv.writer(sgd_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            rgd_out = csv.writer(rgd_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            mgi_out = csv.writer(mgi_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            human_out = csv.writer(human_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')

            # This list is now sorted phylogenetically for the header to be sorted
            out_write_list = [human_out, rgd_out, mgi_out, zfin_out, fb_out, wb_out, sgd_out, sarscov2_out]

            taxon_file_dispatch_dict = {
                'taxid:10116': rgd_out,
                'taxid:9606': human_out,
                'taxid:10090': mgi_out,
                'taxid:6239': wb_out,
                'taxid:559292': sgd_out,
                'taxid:7955': zfin_out,
                'taxid:7227': fb_out,
                'taxid:2697049': sarscov2_out,
                'taxid:4932': sgd_out,
                'taxid:307796': sgd_out,
                'taxid:643680': sgd_out,
                'taxid:574961': sgd_out,
                'taxid:285006': sgd_out,
                'taxid:545124': sgd_out,
                'taxid:764097': sgd_out
            }

            out_to_species_name_dict = {
                rgd_out: 'Rattus norvegicus',
                human_out: 'Homo sapiens',
                mgi_out: 'Mus musculus',
                wb_out: 'Caenorhabditis elegans',
                sgd_out: 'Saccharomyces cerevisiae',
                zfin_out: 'Danio rerio',
                sarscov2_out: 'Severe acute respiratory syndrome coronavirus 2',
                fb_out: 'Drosophila melanogaster'
            }

            out_to_header_taxonid_dict = {
                rgd_out: 'NCBI:txid10116',
                human_out: 'NCBI:txid9606',
                mgi_out: 'NCBI:txid10090',
                wb_out: 'NCBI:txid6239',
                sgd_out: 'NCBI:txid559292',
                zfin_out: 'NCBI:txid7955',
                sarscov2_out: 'NCBI:txid2697049',
                fb_out: 'NCBI:txid7227'
            }

            # Write the comments in the main file.
            filetype = 'Molecular Interactions'
            data_format = 'PSI-MI TAB 2.7 Format'
            database_version = self.context_info.env["ALLIANCE_RELEASE"]
            species_list = []
            taxon_list = []
            for entry in out_write_list:
                taxon_list.append(out_to_header_taxonid_dict[entry])
                species_list.append(out_to_species_name_dict[entry])
            species = ", ".join(species_list)
            taxon_ids = ", ".join(taxon_list)
            taxon_ids = '# TaxonIDs: {}'.format(taxon_ids)
            gen_time = datetime.utcnow().strftime("%Y-%m-%d %H:%M")
            readme = 'https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md'
            response = urllib.request.urlopen(self.context_info.env["HEADER_TEMPLATE_URL"])
            header_template = HeaderTemplate(response.read().decode('ascii'))
            header_dict = {'filetype': filetype, 'data_format': data_format, 'stringency_filter': '',
                           'taxon_ids': taxon_ids, 'database_version': database_version, 'species': species,
                           'gen_time': gen_time, 'readme': readme}
            header = header_template.substitute(header_dict)
            header_rows = [line.strip() for line in header.splitlines() if len(line.strip()) != 0]
            for header_row in header_rows:
                tsvout.writerow([header_row])
            tsvout.writerow(psi_mi_tab_header)

            for entry in out_write_list:
                filetype = 'Molecular Interactions'
                species = out_to_species_name_dict[entry]
                taxon_ids = '# TaxonIDs: {}'.format(out_to_header_taxonid_dict[entry])
                header_dict = {'filetype': filetype, 'data_format': data_format, 'stringency_filter': '',
                               'taxon_ids': taxon_ids, 'database_version': database_version, 'species': species,
                               'gen_time': gen_time, 'readme': readme}
                header = header_template.substitute(header_dict)
                header_rows = [line.strip() for line in header.splitlines() if len(line.strip()) != 0]
                for header_row in header_rows:
                    entry.writerow([header_row])
                entry.writerow(psi_mi_tab_header)

            psi_mi_tab_header.insert(0,'Reason for skipping row.')
            skipped_out.writerow(psi_mi_tab_header)

            # The order of this list is important! Defined in the list above.  Cannot be parallelized
            for filename in parsing_list: 
                logger.info('Parsing file: %s' % (filename))
                filename_type = None
                if filename == imex_filename:
                    filename_type = 'imex'
                elif filename == biogrid_filename:
                    filename_type = 'biogrid'
                elif filename == flybase_filename:
                    filename_type = 'flybase'
                elif filename == wormbase_filename:
                    filename_type = 'wormbase'

                # Declare the tracking dict used to look for duplicates. It tracks sets.
                publication_tracking_dict[filename_type] = set()

                with open(filename, 'r', encoding='utf-8') as tsvin:
                    csv_reader = csv.reader(tsvin, delimiter='\t', quoting=csv.QUOTE_NONE)

                    for row in tqdm(csv_reader):
                        if row[0].startswith("#"):
                            row.insert(0,'Entry starts with # commented out or header')
                            skipped_out.writerow(row)
                            continue

                        if filename_type == 'biogrid':
                            biogrid_interaction_id = re.findall(r'\d+', row[13])[0]
                            # We need to add '-' characters to columns 17-42 for biogrid entries.
                            for _ in range(17,43):
                                row.append('-')
                            # Reassigning values in several columns as described above.
                            row[18] = 'psi-mi:\"MI:0496\"(bait)'
                            row[19] = 'psi-mi:\"MI:0498\"(prey)'
                            row[20] = 'psi-mi:\"MI:0326\"(protein)'
                            row[21] = 'psi-mi:\"MI:0326\"(protein)'
                            if biogrid_interaction_id in self.biogrid_rna_set:
                                row[21] = 'psi-mi:\"MI:0320\"(ribonucleic acid)'
                            row[35] = 'false'

                        try:
                            taxon_id_1 = re.search(r'taxid:\d+', row[9]).group(0)
                        except AttributeError:
                            row.insert(0,'Taxon ID appears to be missing for interactor A from row 10: %s' % row[9])
                            skipped_out.writerow(row)
                            continue # Skip rows where we don't find a taxon entry.

                        try:
                            taxon_id_2 = re.search(r'taxid:\d+', row[10]).group(0)
                        except AttributeError:
                            row.insert(0,'Taxon ID appears to be missing for interactor B from row 11: %s' % row[10])
                            skipped_out.writerow(row)
                            continue # Skip rows where we don't find a taxon entry.

                        if not taxon_id_1 in (taxon_species_set) or not taxon_id_2 in (taxon_species_set):
                            row.insert(0,'Not in Alliance taxon list.')
                            skipped_out.writerow(row)
                            continue # Skip rows where we don't have Alliance species or a blank entry.
                        if taxon_id_1 in possible_yeast_taxon_set: # Change yeast taxon ids to the preferred 'taxid:559292'
                            row[9] = 'taxid:559292(Saccharomyces cerevisiae)'
                        if taxon_id_2 in possible_yeast_taxon_set: # Change yeast taxon ids to the preferred 'taxid:559292'
                            row[10] = 'taxid:559292(Saccharomyces cerevisiae)'

                        row_initial = list(row) # Create a true copy of the list. Used for validation later.

                        # Skip rows with MI:0208 "genetic interaction".
                        if row[11].startswith(interaction_exclusion_set):
                            row.insert(0,'Contains a term from the interaction exclusion set.')
                            skipped_out.writerow(row)
                            continue

                        # Skip rows with undesired interaction types.
                        # Sometimes these columns don't exist in BIOGRID? IndexErrors still write the proper message and skip the entry.
                        if filename_type != 'biogrid': # Biogrid stops at row 16.
                            try:
                                if row[20].startswith(interactor_type_exclusion_set):
                                    row.insert(0,'Contains a term from the interactor type exclusion set.')
                                    skipped_out.writerow(row)
                                    continue
                            except IndexError:
                                    row.insert(0,'Interactor type column not found? Skipping entry.')
                                    skipped_out.writerow(row)
                                    continue
                            try:
                                if row[21].startswith(interactor_type_exclusion_set):
                                    row.insert(0,'Contains a term from the interactor type exclusion set.')
                                    skipped_out.writerow(row)
                                    continue
                            except IndexError:
                                    row.insert(0,'Interactor type column not found? Skipping entry.')
                                    skipped_out.writerow(row)
                                    continue

                            # Skip entries which have 'Expansion method(s)'.
                            if row[15] is not '-':
                                row.insert(0,'Contains an expansion method.')
                                skipped_out.writerow(row)
                                continue


                        interactor_A_resolved, interactor_B_resolved, row = self.resolve_identifiers_by_row(row, mapped_out)

                        if interactor_A_resolved is False and interactor_B_resolved is True:
                            row.insert(0,'Can\'t resolve interactor A identifier, alias, alternate id, or xref against the list of known Alliance identifiers.')
                            skipped_out.writerow(row)
                            continue
                        elif interactor_A_resolved is True and interactor_B_resolved is False:
                            row.insert(0,'Can\'t resolve interactor B identifier, alias, alternate id, or xref against the list of known Alliance identifiers.')
                            skipped_out.writerow(row)
                            continue
                        elif interactor_A_resolved is False and interactor_B_resolved is False:
                            row.insert(0,'Can\'t resolve either interactor A or B identifier, alias, alternate id, or xref against the list of known Alliance identifiers.')
                            skipped_out.writerow(row)
                            continue

                        # Grab the publication information
                        # Also creating a tuple "key" to use for filtering purposes.
                        if row[8] is not None:

                            pubmed_unassigned = re.search(r'pubmed:unassigned', row[8])
                            if pubmed_unassigned is not None:
                                row.insert(0,'Pubmed value is unassigned')
                                skipped_out.writerow(row)
                                continue

                            publication_re = re.search(r'pubmed:\d+', row[8])
                            if publication_re is not None:
                                publication = publication_re.group(0)

                                if publication.startswith('pubmed:88880000'):
                                    biogrid_key = row[13]
                                    publication = self.biogrid_doi_dict[biogrid_key]
                                    row[8] = self.biogrid_doi_dict[biogrid_key]

                                # Capture everything up to the first parenthesis in the taxon column.
                                taxon1 = re.search(r'taxid:\d+', row[9]).group(0)
                                taxon2 = re.search(r'taxid:\d+', row[10]).group(0)
                                # Build a filtering key from the publication, taxon1, and taxon2.
                                tracking_tuple = (publication, taxon1, taxon2)
                                exit_tsv_loop = False
                                for key, value in publication_tracking_dict.items():
                                    if key != filename_type: # Don't look in our current dictionary.
                                        if tracking_tuple in value:
                                            row.insert(0,'Already added this interaction to the export file from %s. Filter criteria: %s' % (key, (tracking_tuple,)))
                                            skipped_out.writerow(row)
                                            exit_tsv_loop = True

                                if exit_tsv_loop == True:
                                    continue

                                # If we loop through all the possible sets and don't continue, add the tuple.
                                publication_tracking_dict[filename_type].add(tracking_tuple)

                        if filename_type != 'biogrid' and len(row) != len(row_initial):
                            print('FATAL ERROR: The row length was changed during processing.')
                            print(row_initial)
                            print(row)
                            quit()

                        # Write the row to all the appropriate files.

                        tsvout.writerow(row)

                        taxon1 = re.search(r'taxid:\d+', row[9]).group(0)
                        taxon2 = re.search(r'taxid:\d+', row[10]).group(0)

                        self.wrote_to_file_already = False


                        if taxon1 == taxon2:
                            try:
                                taxon_file_dispatch_dict[taxon1].writerow(row)
                            except KeyError:
                                pass

                        else:
                            try:
                                taxon_file_dispatch_dict[taxon1].writerow(row)
                            except KeyError:
                                pass

                            try:
                                taxon_file_dispatch_dict[taxon2].writerow(row)
                            except KeyError:
                                pass


    def validate_and_upload_files_to_fms(self):
        logger.info('Summary of files created:')
        logger.info(os.system("ls -alh {}*".format(self.output_dir)))

        upload_location_dict = {
            'alliance_molecular_interactions.tsv': 'COMBINED',
            'alliance_molecular_interactions_fly.tsv': 'FB',
            'alliance_molecular_interactions_sarscov2.tsv': 'SARS-CoV-2',
            'alliance_molecular_interactions_worm.tsv': 'WB',
            'alliance_molecular_interactions_zebrafish.tsv': 'ZFIN',
            'alliance_molecular_interactions_yeast.tsv': 'SGD',
            'alliance_molecular_interactions_rat.tsv': 'RGD',
            'alliance_molecular_interactions_mouse.tsv': 'MGI',
            'alliance_molecular_interactions_human.tsv': 'HUMAN'
        }

        thread_pool = []

        for filename in upload_location_dict.keys():
            upload_location = upload_location_dict[filename]
            compressed_filename = "INTERACTION-MOL_{}.tar.gz".format(upload_location)

            p = multiprocessing.Process(target=self._compress_and_upload, args=(compressed_filename, filename, upload_location))
            p.start()
            thread_pool.append(p)

        Processor.wait_for_threads(thread_pool)


    def _compress_and_upload(self, compressed_filename, filename, upload_location):
        os.chdir(self.output_dir)
        logger.info('Compressing file: {}'.format(filename))
        os.system("tar -czvf {} {}".format(compressed_filename, filename))
        logger.info("tar -czvf {} {}".format(compressed_filename, filename))

        upload_file_prefix = '{}_{}_{}'.format(self.context_info.env['ALLIANCE_RELEASE'], 'INTERACTION-MOL', upload_location)
        file_to_upload = {upload_file_prefix: open(self.output_dir + compressed_filename, 'rb')}

#         self.context_info.env['API_KEY'] = '<insert key here>'      # if don't have have API_KEY in config file, could enter here
        headers = {
            'Authorization': 'Bearer {}'.format(self.context_info.env['API_KEY'])
        }

        logger.info('Attempting upload of data file: {}'.format(compressed_filename))
        logger.debug('Attempting upload with header: {}'.format(headers))
        logger.info("Uploading data to {}) ...".format(self.context_info.env['FMS_API_URL'] + '/api/data/submit/'))

        response = requests.post(self.context_info.env['FMS_API_URL'] + '/api/data/submit/', files=file_to_upload, headers=headers)
        logger.info(response.text)
