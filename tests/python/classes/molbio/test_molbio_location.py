from socialgene.base.molbio import Location


def test_Location_1():
    temp = Location(start=1, end=1, strand=0)
    assert temp.start == 1
    assert temp.end == 1
    assert temp.strand == 0
    assert temp.all_attributes() == {"end": 1, "start": 1, "strand": 0}
    temp = Location(start=100, end=2000, strand=1)
    assert temp.start == 100
    assert temp.end == 2000
    assert temp.strand == 1
    assert temp.all_attributes() == {"end": 2000, "start": 100, "strand": 1}
