import logging, coloredlogs, os, sys, multiprocessing, time
from process import *

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

class PreProcessor(object):

    def pre_process(self):
        logger.info("Running Pre Processor")

        process_dispatch = {
            'MI': MIPreProcess,
            'Interaction': InteractionsPreProcess,
        }

        list_of_process_groups = [
            ['MI', 'Interaction'],
        ]  

        for process_group in list_of_process_groups:
            logger.debug("Process's in group: %s" % process_group)
            thread_pool = []
            for process_name in process_group:
                logger.debug("Process Name: %s" % process_name)
                #config = data_manager.get_config(process_name)
                #logger.debug("Config: %s" % config)
                process = process_dispatch[process_name]()
                p = multiprocessing.Process(target=process.run_pre_process)
                p.start()
                thread_pool.append(p)
            for thread in thread_pool:
                thread.join()

if __name__ == '__main__':
    PreProcessor().pre_process()
