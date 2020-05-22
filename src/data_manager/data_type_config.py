import logging
import sys, os

from transactors import FileTransactor
from common import ContextInfo

from .sub_type_config import SubTypeConfig

logger = logging.getLogger(__name__)


class DataTypeConfig(object):

    def __init__(self, data_type, submission_system_data):
        self.data_type = data_type
        self.submission_system_data = submission_system_data

        self.list_of_subtype_objects = []

    def get_data(self):
        path = 'tmp'
        context_info = ContextInfo()
        if "SAVE_PATH" in context_info.env:
            if context_info.env["SAVE_PATH"]:
                path = context_info.env["SAVE_PATH"]

        # Create our subtype objects.
        for downloadable_item in self.submission_system_data:
            logger.debug("downloadable_item")
            if downloadable_item[2] is not None:
                full_path_to_send = os.path.join(path, downloadable_item[2])
            else:
                full_path_to_send = None  # If we don't have a path.

            sub_type = SubTypeConfig(
                self.data_type, 
                downloadable_item[0], 
                downloadable_item[1], 
                full_path_to_send)

            self.list_of_subtype_objects.append(sub_type)
            logger.debug("data_type_config.py, datatype %s DI0 %s DI1 %s full_path %s end" % (self.data_type, downloadable_item[0], downloadable_item[1], full_path_to_send))

            # Send it off to be queued and executed.
            FileTransactor.execute_transaction(sub_type)

    def check_for_single(self):
        if len(self.list_of_subtype_objects) > 1:
            logger.critical('Called for single item in object containing multiple children.')
            logger.critical('Please check the function calling for this single item.')
            sys.exit(-1)
        else:
            pass

    def get_single_filepath(self):
        self.check_for_single()
        return self.list_of_subtype_objects[0].get_filepath()

    def get_sub_type_objects(self):
        return self.list_of_subtype_objects
