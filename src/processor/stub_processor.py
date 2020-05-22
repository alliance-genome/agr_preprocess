import logging
import multiprocessing

from processor import Processor

logger = logging.getLogger(__name__)


class StubProcessor(Processor):
    def __init__(self, configs):
        super().__init__()
        self.data_type_configs = configs

    def _load_and_process_data(self):
        logger.debug("in StubProcessor")

        source_filepaths = dict()
        stub_config = self.data_type_configs[0]
        for sub_type in stub_config.get_sub_type_objects():
            sub_type_name = sub_type.get_sub_data_type()
            sub_type_filepath = sub_type.get_filepath()
            source_filepaths[sub_type_name] = sub_type_filepath

        for sub_type in source_filepaths:
            logger.debug("Source subtype %s filepath %s" % (sub_type, source_filepaths[sub_type]))

