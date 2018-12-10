import logging, coloredlogs, os, sys, multiprocessing, time

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
        pass


if __name__ == '__main__':
    PreProcessor().pre_process()
