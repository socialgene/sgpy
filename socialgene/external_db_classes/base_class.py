from abc import ABC, abstractmethod
from socialgene.utils.logging import log
import logging
from socialgene.neo4j.neo4j import GraphDriver

logging.getLogger("neo4j").setLevel(logging.WARNING)
logging.getLogger().setLevel(logging.INFO)


class ExternalBaseClass(ABC):
    def _add_to_neo4j(self, statement, database="neo4j", **kwargs):
        summary = (
            GraphDriver()
            .driver.execute_query(
                statement,
                database_=database,
                **kwargs,
            )
            .summary
        )

        if summary.metadata.get("stats"):
            log.info(
                f"{summary.metadata.get('stats').get('properties-set')} properties modified"
            )
        else:
            log.info("No properties modified")
