import csv
import re
import logging
import argparse
import yaml
import os
import time
import coloredlogs
import urllib.request
import json
import glob
import strict_rfc3339
import requests
from retry import retry
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', help='Enable verbose mode.', action='store_true')

args = parser.parse_args()

if args.verbose or ("DEBUG" in os.environ and os.environ['DEBUG'] == "True"):
    debug_level = logging.DEBUG
else:
    debug_level = logging.INFO

coloredlogs.install(level=debug_level,
                    fmt='%(asctime)s %(levelname)s: %(name)s:%(lineno)d: %(message)s',
                    field_styles={
                            'asctime': {'color': 'green'},
                            'hostname': {'color': 'magenta'},
                            'levelname': {'color': 'white', 'bold': True},
                            'name': {'color': 'blue'},
                            'programname': {'color': 'cyan'}
                    })

logger = logging.getLogger(__name__)


class ProcessInteractions(object):

    def __init__(self):
        # List our class variables.
        self.config = None
        self.master_gene_set = set()
        self.master_crossreference_dictionary = dict()
        self.master_crossreference_dictionary['UniProtKB'] = dict()
        self.master_crossreference_dictionary['ENSEMBL'] = dict()
        self.master_crossreference_dictionary['NCBI_Gene'] = dict()
        self.biogrid_ignore_set = set()

        # Load configuration data.
        config_file_loc = 'src/process/config/interactions.yaml'
        logger.debug('Loading configuration file: {}'.format(config_file_loc))
        config_file = open(config_file_loc, 'r')
        self.config = yaml.load(config_file, Loader=yaml.SafeLoader)

        # Look for ENV variables to replace default variables from config file.
        possible_env_keys = ['RELEASE_VERSION', 'FMS_API_URL', 'API_KEY']

        for key in self.config.keys():
            if key in possible_env_keys:
                try:
                    self.config[key] = os.environ[key]
                    logger.info('Found environmental variable for \'{}\'.'.format(key))
                except KeyError:
                    logger.info('Environmental variable not found for \'{}\'. '
                                'Using interactions.yaml value.'.format(key))
                    pass  # If we don't find an ENV variable, keep the value from the config file.

        logger.debug('Initialized with config values: {}'.format(self.config))

    def parse_filepaths_from_s3(self):
        api_url = self.config['FMS_API_URL']
        release_version = self.config['RELEASE_VERSION']

        # Download BGI files from the Alliance.
        for config_entry in self.config['BGI']:
            logger.debug('Obtaining filepath for {}'.format(config_entry))
            query_url = api_url + '/api/datafile/by/{}/{}/{}?latest=true'\
                                  .format(release_version, 'BGI', config_entry)
            with urllib.request.urlopen(query_url) as url:
                data = json.loads(url.read().decode())
                s3_filepath = data[0]['s3Path']
                logger.info('S3 filepath: {}'.format(s3_filepath))
                self.download_file('https://download.alliancegenome.org/' + s3_filepath, config_entry + '.json')

        # Download interaction files from the Alliance, including IMEx and BioGRID.
        list_of_source_files = ['FB', 'WB', 'IMEX', 'BIOGRID', 'BIOGRID-ORGANISM']

        for entry in list_of_source_files:
            logger.info('Obtaining interaction source files from the file management system.')
            query_url = api_url + '/api/datafile/by/{}/{}/{}?latest=true'\
                                  .format(release_version, 'INTERACTION-SOURCE', entry)
            with urllib.request.urlopen(query_url) as url:
                data = json.loads(url.read().decode())
                s3_filepath = data[0]['s3Path']
                logger.info('S3 filepath: {}'.format(s3_filepath))
                self.download_file('https://download.alliancegenome.org/' + s3_filepath, entry + '.zip')
                if entry != 'BIOGRID-ORGANISM':  # This file is unzipped later.
                    logger.info('Extracting file {}.zip with unzip.'.format(entry))
                    os.system('unzip /usr/src/app/download/{}.zip -d /usr/src/app/download/tmp'.format(entry))
                    logger.info('Renaming extracted file.')
                    if entry == 'IMEX':  # Special exception for IMEX because it's 2 files.
                        os.system('mv /usr/src/app/download/tmp/intact.txt /usr/src/app/download/{}.txt'.format(entry))
                        os.system('rm /usr/src/app/download/tmp/*')
                    else:
                        os.system('mv /usr/src/app/download/tmp/* /usr/src/app/download/{}.txt'.format(entry))

    @retry(tries=5, delay=5, logger=logger)
    def download_file(self, query_url, filename):
        directory_and_filename = os.path.join('/usr/src/app/download/', filename)
        logger.info('Attempting to download: {} -> {}'.format(query_url, directory_and_filename))
        urllib.request.urlretrieve(query_url, directory_and_filename)

    def parse_bgi_json(self):
        # We're populating a rather large dictionary to use for looking up Alliance genes by their crossreferences.
        # Edit the list below if you'd like to add more crossreferences to the dictionary.
        # The key of the dictionary is the crossreference and the value is the Alliance gene to which it resolves.
        #
        # We're also populating the "master gene set" for gene lookups later.
        logger.info('Populating master gene set and crossreferences from JSON.')
        for config_entry in self.config['BGI']:
            with open('/usr/src/app/download/' + config_entry + '.json') as json_file:
                data = json.load(json_file)
                logger.info('Scanning {}'.format(config_entry + '.json'))
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

                             # The ids in PSI-MITAB files are lower case, hence the .lower() used above.

        logger.info('Done.')

    def create_rna_protein_ignore(self):
        # Unzip the BioGRID organism file.
        logger.info('Extracting BioGRID organism file with unzip.')
        os.system('unzip /usr/src/app/download/BIOGRID-ORGANISM.zip -d /usr/src/app/download/organism')

        # Generate the ignore list by iterating through the individual organism files.
        for name in glob.glob('/usr/src/app/download/organism/*'):
            for filter_list_item in self.config['BIOGRID_FILTER']:
                if filter_list_item in name:
                    logger.info('Processing filter file: {}'.format(name))
                    with open(name) as tsvfile:
                        reader = csv.reader(tsvfile, delimiter='\t')
                        for row in tqdm(reader):
                            if 'RNA' in row[11]:
                                self.biogrid_ignore_set.add(row[0])


    def resolve_identifiers_by_row(self, row, mapped_out):
        interactor_A_rows = [0, 2, 4, 22]
        interactor_B_rows = [1, 3, 5, 23]

        interactor_A_resolved = False
        interactor_B_resolved = False

        for row_entry in interactor_A_rows:
            try:
                interactor_A_resolved, A_resolved_id = self.resolve_identifier(row[row_entry])
                if interactor_A_resolved is True:
                    break
            except IndexError: # Biogrid has less rows than other files, continue on IndexErrors.
                continue

        for row_entry in interactor_B_rows:
            try:
                interactor_B_resolved, B_resolved_id = self.resolve_identifier(row[row_entry])
                if interactor_B_resolved is True:
                    break
            except IndexError: # Biogrid has less rows than other files, continue on IndexErrors.
                continue

        if A_resolved_id is not None and B_resolved_id is not None:
            mapped_output_rows = [row[13], A_resolved_id, B_resolved_id]
            mapped_out.writerow(mapped_output_rows)

        return interactor_A_resolved, interactor_B_resolved

    def resolve_identifier(self, row_entry):

            list_of_crossref_regex_to_search = [
                'uniprotkb:[\\w\\d_-]*$',
                'ensembl:[\\w\\d_-]*$',
                'entrez gene/locuslink:.*$'
            ]

            # If we're dealing with multiple identifiers separated by a pipe.
            if '|' in row_entry:
                row_entries = row_entry.split('|')
            else:
                row_entries = [row_entry]

            for individual_entry in row_entries:

                # For use in wormbase / flybase lookups.
                # If we run into an IndexError, there's no identifier to resolve and we return False.
                try:
                    entry_stripped = individual_entry.split(':')[1]
                except IndexError:
                    return False, None

                if individual_entry.startswith('uniprotkb:'):
                    individual_entry = individual_entry.split('-')[0]

                prefixed_identifier = None

                if entry_stripped.startswith('WB'): # TODO implement regex for WB / FB gene identifiers.
                    prefixed_identifier = 'WB:' + entry_stripped
                    if prefixed_identifier in self.master_gene_set:
                        return True, prefixed_identifier
                    else:
                        return False, None
                elif entry_stripped.startswith('FB'): # TODO implement regex for WB / FB gene identifiers.
                    prefixed_identifier = 'FB:' + entry_stripped
                    if prefixed_identifier in self.master_gene_set:
                        return True, prefixed_identifier
                    else:
                        return False, None

                for regex_entry in list_of_crossref_regex_to_search:
                    regex_output = re.findall(regex_entry, individual_entry)
                    if regex_output is not None:
                        for regex_match in regex_output: # We might have multiple regex matches. Search them all against our crossreferences.
                            identifier = regex_match
                            for crossreference_type in self.master_crossreference_dictionary.keys():
                                # Using lowercase in the identifier to be consistent with Alliance lowercase identifiers.
                                if identifier.lower() in self.master_crossreference_dictionary[crossreference_type]:
                                    return True, identifier.lower() # Return 'True' if we find an entry.

            # If we can't resolve any of the crossReferences, return None
            return False, None

    def get_data(self):

        imex_filename = '/usr/src/app/download/IMEX.txt'
        wormbase_filename = '/usr/src/app/download/WB.txt'
        flybase_filename = '/usr/src/app/download/FB.txt'
        biogrid_filename = '/usr/src/app/download/BIOGRID.txt'

        # The order of this list is important.
        parsing_list = [wormbase_filename, flybase_filename, biogrid_filename, imex_filename]

        # TODO Taxon species needs to be pulled out into a standalone module to be used by other scripts.
        # TODO External configuration script for these types of filters? Not a fan of hard-coding.
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
        with open('output/alliance_molecular_interactions.tsv', 'w', encoding='utf-8') as tsvout, \
             open('output/skipped_entries.txt', 'w', encoding='utf-8') as skipped_out, \
             open('output/alliance_molecular_interactions_fly.tsv', 'w', encoding='utf-8') as fb_out, \
             open('output/alliance_molecular_interactions_worm.tsv', 'w', encoding='utf-8') as wb_out, \
             open('output/alliance_molecular_interactions_zebrafish.tsv', 'w', encoding='utf-8') as zfin_out, \
             open('output/alliance_molecular_interactions_yeast.tsv', 'w', encoding='utf-8') as sgd_out, \
             open('output/alliance_molecular_interactions_rat.tsv', 'w', encoding='utf-8') as rgd_out, \
             open('output/alliance_molecular_interactions_mouse.tsv', 'w', encoding='utf-8') as mgi_out, \
             open('output/alliance_molecular_interactions_human.tsv', 'w', encoding='utf-8') as human_out, \
             open('output/mapped_entries.txt', 'a+', encoding='utf-8') as mapped_out:

            mapped_out = csv.writer(mapped_out, quotechar = '', quoting=csv.QUOTE_NONE, delimiter='\t')
            tsvout = csv.writer(tsvout, quotechar = '', quoting=csv.QUOTE_NONE, delimiter='\t')
            skipped_out = csv.writer(skipped_out, quotechar = '', quoting=csv.QUOTE_NONE, delimiter='\t')
            fb_out = csv.writer(fb_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            wb_out = csv.writer(wb_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            zfin_out = csv.writer(zfin_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            sgd_out = csv.writer(sgd_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            rgd_out = csv.writer(rgd_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            mgi_out = csv.writer(mgi_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')
            human_out = csv.writer(human_out, quotechar='', quoting=csv.QUOTE_NONE, delimiter='\t')

            out_write_list = [fb_out, wb_out, zfin_out, sgd_out, rgd_out, mgi_out, human_out]

            taxon_file_dispatch_dict = {
                'taxid:10116': rgd_out,
                'taxid:9606': human_out,
                'taxid:10090': mgi_out,
                'taxid:6239': wb_out,
                'taxid:559292': sgd_out,
                'taxid:7955': zfin_out,
                'taxid:7227': fb_out,
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
                fb_out: 'Drosophila melanogaster'
            }

            # Write the comments in the main file.
            the_time = strict_rfc3339.now_to_rfc3339_localoffset()
            tsvout.writerow(['# Molecular Interactions'])
            tsvout.writerow(['# Alliance of Genome Resources'])
            tsvout.writerow(['# alliancegenome.org'])
            tsvout.writerow(['# alliance-helpdesk@lists.stanford.edu'])
            tsvout.writerow(['#'])
            tsvout.writerow(['# PSI-MI TAB 2.7 Format'])
            tsvout.writerow(['# Alliance Version {}'.format(self.config['RELEASE_VERSION'])])
            tsvout.writerow(['# Date Produced: {}'.format(the_time)])
            tsvout.writerow(['#'])
            # Write the headers
            tsvout.writerow(psi_mi_tab_header)

            for entry in out_write_list:
                the_time = strict_rfc3339.now_to_rfc3339_localoffset()
                entry.writerow(['# Molecular Interactions for {}'.format(out_to_species_name_dict[entry])])
                entry.writerow(['# Alliance of Genome Resources'])
                entry.writerow(['# alliancegenome.org'])
                entry.writerow(['# alliance-helpdesk@lists.stanford.edu'])
                entry.writerow(['#'])
                entry.writerow(['# PSI-MI TAB 2.7 Format'])
                entry.writerow(['# Alliance Version {}'.format(self.config['RELEASE_VERSION'])])
                entry.writerow(['# Date Produced: {}'.format(the_time)])
                entry.writerow(['#'])
                # Write the headers
                entry.writerow(psi_mi_tab_header)

            psi_mi_tab_header.insert(0,'Reason for skipping row.')
            skipped_out.writerow(psi_mi_tab_header)

            for filename in parsing_list: # The order of this list is important! Defined in the list above.
                print('Parsing file: %s' % (filename))
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
                    tsvin = csv.reader(tsvin, delimiter='\t', quoting=csv.QUOTE_NONE)
                    if filename_type is not 'wormbase':
                        next(tsvin, None) # Skip the headers (no headers in wormbase file).

                    for row in tqdm(tsvin):

                        # Skip entries with unassigned pubmed IDs.
                        if 'pubmed:unassigned' in row[8]:
                            row.insert(0,'Entry contains a pubmed:unassigned identifier.')
                            skipped_out.writerow(row)
                            continue

                        if filename_type == 'biogrid':
                            biogrid_interaction_id = re.findall(r'\d+', row[13])[0]
                            # Exclude interations from the biogrid_prot_rna_set.
                            if biogrid_interaction_id in self.biogrid_ignore_set:
                                row.insert(0,'Entry appears in biogrid protein-rna filter set.')
                                skipped_out.writerow(row)
                                continue
                            # We need to add '-' characters to columns 16-42 for biogrid entries.
                            for _ in range(16,43):
                                row.append('-')
                            # Reassigning values in several columns as described above.
                            row[18] = 'psi-mi:\"MI:0496\"(bait)'
                            row[19] = 'psi-mi:\"MI:0498\"(prey)'
                            row[20] = 'psi-mi:\"MI:0326\"(protein)'
                            row[21] = 'psi-mi:\"MI:0326\"(protein)'

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

                        interactor_A_resolved, interactor_B_resolved = self.resolve_identifiers_by_row(row, mapped_out)

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
                            publication_re = re.search(r'pubmed:\d+', row[8])
                            if publication_re is not None:
                                publication = publication_re.group(0)
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

                        try:
                            taxon_file_dispatch_dict[taxon1].writerow(row)
                            self.wrote_to_file_already = True
                        except KeyError:
                            pass

                        try:
                            if self.wrote_to_file_already is False:
                                taxon_file_dispatch_dict[taxon2].writerow(row)
                        except KeyError:
                            pass

    def validate_and_upload_files_to_fms(self):

        logger.info('Summary of files created:')
        logger.info(os.system('ls -alh /usr/src/app/output/*'))

        upload_location_dict = {
            'alliance_molecular_interactions.tsv': 'COMBINED',
            'alliance_molecular_interactions_fly.tsv': 'FB',
            'alliance_molecular_interactions_worm.tsv': 'WB',
            'alliance_molecular_interactions_zebrafish.tsv': 'ZFIN',
            'alliance_molecular_interactions_yeast.tsv': 'SGD',
            'alliance_molecular_interactions_rat.tsv': 'RGD',
            'alliance_molecular_interactions_mouse.tsv': 'MGI',
            'alliance_molecular_interactions_human.tsv': 'HUMAN'
        }

        # time.sleep(100000)  # Stop the container from exiting so we can check the files.

        for filename in upload_location_dict.keys():

            logger.info('Compressing file: {}'.format(filename))
            compressed_filename = "{}.tar.gz".format(upload_location_dict[filename])

            os.chdir('/usr/src/app/output/')
            os.system("tar -czvf {} {}".format(compressed_filename, filename))
            logger.info(os.system('ls -alh /usr/src/app/output/*'))

            upload_file_prefix = '{}_{}_{}'.\
                format(self.config['RELEASE_VERSION'], 'INTERACTION-MOL', upload_location_dict[filename])

            file_to_upload = {upload_file_prefix: open('/usr/src/app/output/' + compressed_filename, 'rb')}

            headers = {
                'Authorization': 'Bearer {}'.format(self.config['API_KEY'])
            }

            logger.info('Attempting upload of data file: {}'.format(compressed_filename))
            logger.debug('Attempting upload with header: {}'.format(headers))
            logger.info("Uploading data to {}) ...".format(self.config['FMS_API_URL'] + '/api/data/submit/'))

            response = requests.post(self.config['FMS_API_URL'] + '/api/data/submit/', files=file_to_upload,
                                     headers=headers)
            logger.info(response.text)


if __name__ == "__main__":
    interactions = ProcessInteractions()
    interactions.parse_filepaths_from_s3()
    interactions.parse_bgi_json()
    interactions.create_rna_protein_ignore()
    interactions.get_data()
    interactions.validate_and_upload_files_to_fms()
