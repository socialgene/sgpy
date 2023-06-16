"""
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.description) = "True" SET p1:description;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.partial_on_complete_genome) = "True" SET p1:partial_on_complete_genome;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.missing_start) = "True" SET p1:missing_start;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.missing_stop) = "True" SET p1:missing_stop;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.internal_stop) = "True" SET p1:internal_stop;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.partial_in_the_middle_of_a_contig) = "True" SET p1:partial_in_the_middle_of_a_contig;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.missing_N_terminus) = "True" SET p1:missing_N_terminus;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.missing_C_terminus) = "True" SET p1:missing_C_terminus;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.frameshifted) = "True" SET p1:frameshifted;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.too_short_partial_abutting_assembly_gap) = "True" SET p1:too_short_partial_abutting_assembly_gap;
MATCH n=()-[e:ENCODES]->(p1:protein) WHERE (e.incomplete) = "True" SET p1:incomplete;
"""
