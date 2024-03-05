MATCH (n:assembly) where  n.lat_lon is not null
WITH n, split(n.lat_lon, ' ') AS parts
WITH n, toFloat(parts[0]) AS latitude, parts[1] AS lat_direction, toFloat(parts[2]) AS longitude, parts[3] AS lon_direction
WITH n, 
  CASE lat_direction WHEN 'N' THEN latitude WHEN 'S' THEN -latitude ELSE null END AS latitude_value,
  CASE lon_direction WHEN 'E' THEN longitude WHEN 'W' THEN -longitude ELSE null END AS longitude_value
SET n.geolocation= point({latitude: latitude_value, longitude: longitude_value})
//return   point({latitude: latitude_value, longitude: longitude_value}) AS p limit 2
//RETURN n.latitude, n.longitude limit 2





MATCH z1=(n:pfam {name:"Trp_halogenase"})<-[:SOURCE_DB]-(h1:hmm)
MATCH z2=(h1)-[:ANNOTATES]-(:protein)<-[e1:ENCODES]-(n1:nucleotide)
CALL {
    WITH n1, e1
    MATCH z3=(an1:antismash)<-[:SOURCE_DB]-(:hmm)-[:ANNOTATES]->(p1:protein)<-[e2:ENCODES]-(n1)
    MATCH z4=(an2:antismash)<-[:SOURCE_DB]-(:hmm)-[:ANNOTATES]->(p1)
    WHERE an1.name ="Condensation"
        AND an2.name in ["AMP-binding", "A-OX"] AND abs(e1.start - e2.start) < 10000 AND e1.strand = e2.strand
    MATCH z5=(:amrfinder)<-[:SOURCE_DB]-(:hmm)-[:ANNOTATES]->(p2:protein)<-[e3:ENCODES]-(n1)
    WHERE abs(e1.start - e3.start) < 50000 AND e1.strand = e3.strand
    MATCH (n1)-[:ASSEMBLES_TO]->(az:assembly)
    WHERE az.lat_lon is not null
    return az
 } in transactions of 1 rows
 RETURN distinct az
