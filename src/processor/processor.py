import logging
import os, sys
import time
from common import ContextInfo
import requests
import gzip
import shutil

logger = logging.getLogger(__name__)


class Processor(object):

    def __init__(self):
        self.context_info = ContextInfo()

    def run_processor(self):
        self._load_and_process_data()

    def wait_for_threads(thread_pool, queue=None):
        logger.debug("Waiting for Threads to finish: %s" % len(thread_pool))

        while len(thread_pool) > 0:
            logger.debug("Checking Threads: %s" % len(thread_pool))
            for (index, thread) in enumerate(thread_pool):
                logger.debug("Thread Alive: %s Exitcode: %s" % (thread.is_alive(), thread.exitcode))
                if (thread.exitcode is None or thread.exitcode == 0) and not thread.is_alive():
                    logger.debug("Thread Finished Removing from pool: ")
                    thread.join()
                    del thread_pool[index]
                elif thread.exitcode is not None and thread.exitcode != 0:
                    logger.debug("Thread has Problems Killing Children: ")
                    for thread1 in thread_pool:
                        thread1.terminate()
                    sys.exit(-1)
                else:
                    pass

            if queue is not None:
                logger.debug("Queue Size: %s" % queue.qsize())
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
