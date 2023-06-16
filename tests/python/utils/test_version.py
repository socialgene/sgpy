from importlib.metadata import version

from socialgene.utils.version import main


def test_version_print(capsys):
    main()
    captured = capsys.readouterr()
    # print adds newline so strip()
    assert captured.out.strip() == f"{version('socialgene')}"
