import pytest

from socialgene.hashing.hashing import hash_aminos, hasher, sha512t24u, use_hashlib


@pytest.mark.parametrize(
    "n,e",
    [
        ("a", "H0D8ktokFpR1CXnubPWC8tXX0o4YM13g"),
        ("bb", "JLGoEtTjU1wGARxDCquj9Z0y82Jj3cuZ"),
        (" ", "-Q3dd-QA3-aj_PR5sAse4p5wFcW7jNcP"),
        ("\n", "vmiIOMqGhuXJBom_KrWFzvETfJmbSMcL"),
        (str("\\sdsdssdsd\335"), "jvIIYScQeglIZSLg2dAgUwvWwTr9V5xq"),  # noqa: W605
    ],
)
def test_sha512t24u(n, e):
    assert sha512t24u(n) == e


def test_sha512t24_fail():
    with pytest.raises(Exception):
        sha512t24u(1)
    with pytest.raises(Exception):
        sha512t24u([])
    with pytest.raises(Exception):
        sha512t24u({})
    with pytest.raises(Exception):
        sha512t24u(())
    with pytest.raises(Exception):
        sha512t24u(["hi"])


def test_hash_aminos():
    # hash_aminos just makes sure input is all uppercase
    assert sha512t24u("A") == hash_aminos("A")
    assert sha512t24u("A") == hash_aminos("a")
    assert hash_aminos("MTTQAPTFTQ") == "gKM7C-cUyDhPyIFv1xwB-zZzUIRLATmc"


@pytest.mark.parametrize(
    "n,e",
    [
        (
            "sha512",
            "1f40fc92da241694750979ee6cf582f2d5d7d28e18335de05abc54d0560e0f5302860c652bf08d560252aa5e74210546f369fbbbce8c12cfc7957b2652fe9a75",
        ),
        ("sha256", "ca978112ca1bbdcafac231b39a23dc4da786eff8147c4e72b9807785afee48bb"),
        (
            "sha384",
            "54a59b9f22b0b80880d8427e548b7c23abd873486e1f035dce9cd697e85175033caa88e6d57bc35efae0b5afd3145f31",
        ),
        ("md5", "0cc175b9c0f1b6a831c399e269772661"),
        ("sha224", "abd37534c7d9a2efb9465de931cd7055ffdb8879563ae98078d6d6d5"),
    ],
)
def test_use_hashlib(n, e):
    for i in ("sha512", "sha256", "sha384", "md5", "sha224"):
        assert use_hashlib("a", n) == e


@pytest.mark.parametrize(
    "n,e",
    [
        (
            "sha512",
            "1f40fc92da241694750979ee6cf582f2d5d7d28e18335de05abc54d0560e0f5302860c652bf08d560252aa5e74210546f369fbbbce8c12cfc7957b2652fe9a75",
        ),
        ("sha256", "ca978112ca1bbdcafac231b39a23dc4da786eff8147c4e72b9807785afee48bb"),
        (
            "sha384",
            "54a59b9f22b0b80880d8427e548b7c23abd873486e1f035dce9cd697e85175033caa88e6d57bc35efae0b5afd3145f31",
        ),
        ("md5", "0cc175b9c0f1b6a831c399e269772661"),
        ("sha224", "abd37534c7d9a2efb9465de931cd7055ffdb8879563ae98078d6d6d5"),
    ],
)
def test_hasher(n, e):
    assert hasher("a", n) == e
    assert hasher("a", "seguid") == "bc1M4j2I4u6VaLpUbAB8Y9kTHBs"
