import logging, coloredlogs, os, multiprocessing, time, argparse, time

from processor import *
from common import ContextInfo  # Must be the last timeport othersize program fails
from data_manager import DataFileManager
from transactors import FileTransactor



parser = argparse.ArgumentParser(description='Load data into the Neo4j database for the Alliance of Genome Resources.')
parser.add_argument('-c', '--config', help='Specify the filename of the YAML config. It must reside in the src/config/ directory', default='default.yml')
parser.add_argument('-v', '--verbose', help='Enable DEBUG mode for logging.', action='store_true')
args = parser.parse_args()

# set context info
context_info = ContextInfo()
context_info.config_file_location = os.path.abspath('src/config/' + args.config)
if args.verbose:
    context_info.env["DEBUG"] = True

debug_level = logging.DEBUG if context_info.env["DEBUG"] else logging.INFO

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
logging.getLogger("preprocessor").setLevel(logging.ERROR)


class AggregatePreprocessor(object):

    def run_preprocessor(self):

        if args.verbose:
            logger.warn('DEBUG mode enabled!')
            time.sleep(3)

        start_time = time.time()
        data_manager = DataFileManager(context_info.config_file_location)
        logger.info("config_file_location %s" % (context_info.config_file_location))

        ft = FileTransactor()

        ft.start_threads(data_manager.get_FT_thread_settings())
        data_manager.download_and_validate()
        logger.info("finished downloading now doing thread cleanup")
        ft.check_for_thread_errors()
        logger.info("finished threads waiting for queues")
        ft.wait_for_queues()

        logger.info("finished queues waiting for shutdown")
        ft.shutdown()

        configs_dict = {
            'INTERACTION-SOURCE-MOL': ['INTERACTION-SOURCE','BGI'],
            'INTERACTION-SOURCE-GEN': ['INTERACTION-SOURCE','BGI']
        }

        config_dict = {
            'INTERACTION-SOURCE-MOL': 'INTERACTION-SOURCE',
            'INTERACTION-SOURCE-GEN': 'INTERACTION-SOURCE'
        }


        processor_dispatch = {
            'INTERACTION-SOURCE-MOL': InteractionMolecularProcessor,
            'INTERACTION-SOURCE-GEN': InteractionGeneticProcessor
        }

        list_of_processor_groups = [
            ['INTERACTION-SOURCE-MOL','INTERACTION-SOURCE-GEN']
        ]

        processor_time_tracker_list = []

        for processor_group in list_of_processor_groups:
            processor_group_start_time = time.time()
            logger.info("Starting Processor group: %s" % processor_group)
            thread_pool = []
            for processor_name in processor_group:
                logger.info("Processor Name: %s" % processor_name)

                configs = []
                for config_type in configs_dict[processor_name]:
                    config = data_manager.get_config(config_type)
                    if config is not None:
                        configs.append(config)
                    else:
                        logger.info("No Config found for: %s %s" % (processor_name, config_type))

                if len(configs) > 0:
                    processor = processor_dispatch[processor_name](configs)
                    p = multiprocessing.Process(target=processor.run_processor)
                    p.start()
                    thread_pool.append(p)
                else:
                    logger.info("No Configs found for: %s" % processor_name)

            Processor.wait_for_threads(thread_pool)

            logger.info("Waiting for Queues to sync up")
            processor_elapsed_time = time.time() - processor_group_start_time
            processor_time_message = ("Finished Processor group: %s, Elapsed time: %s" % (processor_group, time.strftime("%H:%M:%S", time.gmtime(processor_elapsed_time))))
            logger.info(processor_time_message)
            processor_time_tracker_list.append(processor_time_message)

        end_time = time.time()
        elapsed_time = end_time - start_time

        for time_item in processor_time_tracker_list:
            logger.info(time_item)
        logger.info('PreProcess finished. Elapsed time: %s' % time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))


if __name__ == '__main__':
    AggregatePreprocessor().run_preprocessor()

