import logging, yaml, os, sys, json, urllib3, requests

from cerberus import Validator

from files import JSONFile
from common import Singleton
from common import ContextInfo
from .data_type_config import DataTypeConfig

logger = logging.getLogger(__name__)


class DataFileManager(metaclass=Singleton):
    
    def __init__(self, config_file_loc):

        self.context_info = ContextInfo()

        # Load config yaml.
        logger.debug('Loading config file: %s' % config_file_loc)
        config_file = open(config_file_loc, 'r')
        self.config_data = yaml.load(config_file, Loader=yaml.SafeLoader)
        logger.debug("Config Data: %s" % self.config_data)

        # Load validation yaml.
        validation_yaml_file_loc = os.path.abspath('src/config/validation.yml')
        logger.debug('Loading validation schema: %s' % validation_yaml_file_loc)
        validation_schema_file = open(validation_yaml_file_loc, 'r')
        self.validation_schema = yaml.load(validation_schema_file, Loader=yaml.SafeLoader)

        # Assign values for thread counts.
        self.FileTransactorThreads = self.config_data['FileTransactorThreads']

        # Loading a JSON blurb from a file as a placeholder for submission system query.
        other_file_meta_data = os.path.abspath('src/config/local_submission.json')
        self.non_submission_system_data = JSONFile().get_data(other_file_meta_data)
        urllib3.disable_warnings()
        self.http = urllib3.PoolManager()

        # use the recently created snapshot
        api_url = self.context_info.env["FMS_API_URL"] + '/api/snapshot/release/' + self.context_info.env["ALLIANCE_RELEASE"]
        logger.info(api_url)

        submission_data = self.http.request('GET', api_url)

        if submission_data.status != 200:
            logger.error("Status: %s" % submission_data.status)
            logger.error("No Data came from API: %s" % api_url)
            sys.exit(-1)

        self.snapshot_submission_system_data = json.loads(submission_data.data.decode('UTF-8'))
        logger.debug(self.snapshot_submission_system_data)

        for dataFile in self.non_submission_system_data['snapShot']['dataFiles']:
            self.snapshot_submission_system_data['snapShot']['dataFiles'].append(dataFile)

        logger.debug(self.snapshot_submission_system_data)

        # List used for MOD and data type objects.
        self.master_data_dictionary = {}

        # Dictionary for transformed submission system data.
        self.transformed_submission_system_data = {}

        # process config file during initialization
        self.process_config()
        
    def get_FT_thread_settings(self):
        return self.FileTransactorThreads

    def get_config(self, data_type):
        # Get the object for a data type. If the object doesn't exist, this returns None.
        logger.debug("Getting config for: [%s] -> Config[%s]" % (data_type, self.master_data_dictionary))
        return self.master_data_dictionary.get(data_type)

    def dispatch_to_object(self):
        # This function sends off our data types to become DataTypeConfig objects.
        # The smaller SubTypeConfig objects are created in the DataTypeConfig functions, see data_type_config.py.
        
        for config_entry in self.transformed_submission_system_data.keys():
            # Skip string entries (e.g. schemaVersion, releaseVersion).
            if isinstance(self.transformed_submission_system_data[config_entry], str):
                continue

            logger.debug('Processing DataType: %s' % config_entry)

            # Create our data type object and add it to our master dictionary filed under the config_entry.
            # e.g. Create BGI DataTypeConfig object and file it under BGI in the dictionary.
            self.master_data_dictionary[config_entry] = DataTypeConfig(config_entry,
                                                                       self.transformed_submission_system_data[config_entry])

    def download_and_validate(self):
        logger.debug('Beginning download and validation.')
        for entry in self.master_data_dictionary.keys():
            logger.debug('Downloading %s data.' % entry)
            if isinstance(self.master_data_dictionary[entry], DataTypeConfig):  # If we're dealing with an object.
                self.master_data_dictionary[entry].get_data()
                logger.debug('done with %s data.' % entry)

    def process_config(self):
        # This checks for the validity of the YAML file.
        # See src/config/validation.yml for the layout of the schema.
        # TODO Add requirement checking and more validation to the YAML schema.

        validator = Validator(self.validation_schema)
        validation_results = validator.validate(self.config_data)

        if validation_results is True:
            logger.debug('Config file validation successful.')
        else:
            logger.critical('Config file validation unsuccessful!')
            for field, values in validator.errors.items():
                for value in values:  # May have more than one error per field.
                    message = field + ': ' + value
                    logger.critical(message)
            logger.critical('Exiting')
            sys.exit(-1)

        # Query the submission system for the required data.
        self.query_submission_system()

        # Create our DataTypeConfig (which in turn create our SubTypeConfig) objects.
        self.dispatch_to_object()

    def _search_submission_data(self, dataType, dataSubType):

            try:
                returned_dict = next(item for item in self.snapshot_submission_system_data['snapShot']['dataFiles']
                                     if item['dataType'].get('name') == dataType and item['dataSubType'].get('name') == dataSubType)

            except StopIteration:
                logger.debug('dataType: %s subType: %s not found in submission system data.' % (dataType, dataSubType))
                logger.debug('Creating entry with \'None\' path and extracted path.')
                returned_dict = {
                    'dataType': dataType,
                    'subType': dataSubType,
                    'path': None,
                    'tempExtractedFile': None
                }

            return returned_dict

    def _query_api_datafile_latest(self, dataType, dataSubType):
        api_url = self.context_info.env["FMS_API_URL"] + '/api/datafile/by/' + self.context_info.env["ALLIANCE_RELEASE"] + '/' + dataType + '/' + dataSubType + '?latest=true'
        logger.debug(api_url)

        submission_data = self.http.request('GET', api_url)
        if submission_data.status != 200:
            logger.error("Status: %s" % submission_data.status)
            logger.error("No Data came from API: %s" % api_url)
            sys.exit(-1)

        endpoint_submission_system_data = json.loads(submission_data.data.decode('UTF-8'))
        logger.debug(endpoint_submission_system_data)

        s3Path = endpoint_submission_system_data[0].get('s3Path')
        returned_dict = {
            'dataType': dataType,
            'subType': dataSubType,
            's3Path': s3Path,
            'tempExtractedFile': None
        }
        return returned_dict


    def query_submission_system(self):

        self.transformed_submission_system_data['releaseVersion'] = self.snapshot_submission_system_data['snapShot']['releaseVersion']['releaseVersion']

        config_values_to_ignore = [
            'releaseVersion',  # Manually assigned above.
            'schemaVersion',  # There is no endpoint for latest schema version in api
            'FileTransactorThreads',
            'Neo4jTransactorThreads'
        ]

        for datatype in self.config_data.keys():  # Iterate through our config file.
            logger.debug("Datatype: %s" % datatype)
            if datatype not in config_values_to_ignore:  # Skip these entries.
                self.transformed_submission_system_data[datatype] = []  # Create our empty list.
                for sub_datatype in self.config_data[datatype]:
                      # to process by querying the api for the latest path 
                    submission_system_dict = self._query_api_datafile_latest(datatype, sub_datatype)
                      # to process by using the release snapshot for that path
#                   submission_system_dict = self._search_submission_data(datatype, sub_datatype)

                    path = submission_system_dict.get('s3Path')
                    logger.debug("datatype %s sub_datatype %s path %s" % (datatype, sub_datatype, path))
                    tempExtractedFile = submission_system_dict.get('tempExtractedFile')
                    logger.debug("tempExtractedFile %s" % tempExtractedFile)
                    if tempExtractedFile is None or tempExtractedFile == '':
                        tempExtractedFile = submission_system_dict.get('s3Path')

                    self.transformed_submission_system_data[datatype].append([sub_datatype, path, tempExtractedFile])
            else:
                logger.debug("Ignoring datatype: %s" % datatype)
                        
        logger.debug("Loaded Types: %s" % self.transformed_submission_system_data)
