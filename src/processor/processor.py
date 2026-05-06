import logging
import os, sys
import time
from common import ContextInfo
import requests
import gzip
import shutil
import boto3
from botocore.credentials import InstanceMetadataProvider
from botocore.utils import InstanceMetadataFetcher
from botocore.session import Session as BotocoreSession

logger = logging.getLogger(__name__)


class Processor(object):

    def __init__(self):
        self.context_info = ContextInfo()

    def run_processor(self):
        self._load_and_process_data()

    def wait_for_threads(thread_pool, queue=None):
        logger.info("Waiting for Threads to finish: %s" % len(thread_pool))

        while len(thread_pool) > 0:
            logger.info("Checking Threads: %s" % len(thread_pool))
            for (index, thread) in enumerate(thread_pool):
                logger.info("Thread Alive: %s Exitcode: %s" % (thread.is_alive(), thread.exitcode))
                if (thread.exitcode is None or thread.exitcode == 0) and not thread.is_alive():
                    logger.info("Thread Finished Removing from pool: ")
                    thread.join()
                    del thread_pool[index]
                elif thread.exitcode is not None and thread.exitcode != 0:
                    logger.info("Thread has Problems Killing Children: ")
                    for thread1 in thread_pool:
                        thread1.terminate()
                    sys.exit(-1)
                else:
                    pass

            if queue is not None:
                logger.info("Queue Size: %s" % queue.qsize())
                if queue.empty():
                    queue.join()
                    return
            time.sleep(5)

    def process_query_params(self, query_list_with_params):
        # generators = list of yielded lists from parser
        # query_list_with_parms = list of queries, each with batch size and CSV file name.
        query_and_file_names = []

        for query_params in query_list_with_params:
            cypher_query_template = query_params.pop(0)  # Remove the first query + batch size + CSV file name
            #  from the list. Format the query with all remaining paramenters.
            query_to_run = cypher_query_template % tuple(query_params)

            while len(query_params) > 2:  # We need to remove extra params before we append
                # the modified query. Assuming the last entry in the list is the filepath
                query_params.pop()

            file_name = query_params.pop()
            query_and_file_names.append([query_to_run, file_name])

        return query_and_file_names

    def fms_upload(self, dataType, dataSubType, filepath_uncompressed):

        filepath_compressed = filepath_uncompressed + ".gz"

        with open(self.output_dir + filepath_uncompressed, 'rb') as f_in:
            with gzip.open(self.output_dir + filepath_compressed, 'wb') as f_out:
               shutil.copyfileobj(f_in, f_out)

        upload_file_prefix = '{}_{}_{}'.format(self.context_info.env['ALLIANCE_RELEASE'], dataType, dataSubType)

        file_to_upload = {upload_file_prefix: open(self.output_dir + filepath_compressed, 'rb')}

#       self.context_info.env['API_KEY'] = '<insert key here>'      # if don't have have API_KEY in config file, could enter here
        headers = {
            'Authorization': 'Bearer {}'.format(self.context_info.env['API_KEY'])
        }

        logger.info('Attempting upload of data file: {}'.format(filepath_compressed))
        logger.info('Attempting upload with header: {}'.format(headers))
        logger.info("Uploading data to %s %s %s) ...", upload_file_prefix, filepath_uncompressed, self.context_info.env['FMS_API_URL'] + '/api/data/submit/')

        response = requests.post(self.context_info.env['FMS_API_URL'] + '/api/data/submit/', files=file_to_upload, headers=headers)
        logger.info(response.text)

        self.s3_upload(dataType, dataSubType, filepath_uncompressed, filepath_compressed)

    def s3_upload(self, dataType, dataSubType, filepath_uncompressed, filepath_compressed):
        bucket = self.context_info.env['S3_BUCKET']
        release = self.context_info.env['ALLIANCE_RELEASE']
        local_path = self.output_dir + filepath_compressed

        ext = os.path.splitext(filepath_uncompressed)[1].lstrip('.').lower()
        format_label = ext.upper()
        s3_filename = '{}_{}_{}.{}.gz'.format(dataType, format_label, dataSubType, ext)
        s3_key = '{}/downloads/{}'.format(release, s3_filename)

        s3_client = self._build_s3_client()

        logger.info('Uploading %s to s3://%s/%s', local_path, bucket, s3_key)
        s3_client.upload_file(local_path, bucket, s3_key)
        logger.info('Uploaded s3://%s/%s', bucket, s3_key)

    def _build_s3_client(self):
        profile = self.context_info.env.get('AWS_PROFILE')
        access_key = self.context_info.env.get('AWS_ACCESS_KEY')
        secret_key = self.context_info.env.get('AWS_SECRET_KEY')

        if profile:
            logger.info('Using AWS profile: %s', profile)
            return boto3.Session(profile_name=profile).client('s3')

        if access_key and secret_key:
            logger.info('Using static AWS access keys')
            return boto3.client('s3', aws_access_key_id=access_key, aws_secret_access_key=secret_key)

        instance_creds = InstanceMetadataProvider(iam_role_fetcher=InstanceMetadataFetcher(timeout=2, num_attempts=2)).load()
        if instance_creds is not None:
            logger.info('Using EC2 instance profile credentials')
            botocore_session = BotocoreSession()
            botocore_session._credentials = instance_creds
            return boto3.Session(botocore_session=botocore_session).client('s3')

        logger.info('Using default AWS credential chain')
        return boto3.client('s3')
