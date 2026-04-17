from __future__ import annotations

import types

import pytest

from gatk_sv_ploidy import cli


def test_print_usage_mentions_subcommands(capsys) -> None:
    cli._print_usage()
    captured = capsys.readouterr()
    assert "Whole-genome aneuploidy detection" in captured.out
    assert "preprocess" in captured.out
    assert "infer" in captured.out
    assert "pull-snps" in captured.out


def test_main_exits_with_usage_when_no_subcommand(monkeypatch, capsys) -> None:
    monkeypatch.setattr(cli.sys, "argv", ["gatk-sv-ploidy"])

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 1
    assert "Usage:" in capsys.readouterr().out


def test_main_exits_zero_for_help(monkeypatch, capsys) -> None:
    monkeypatch.setattr(cli.sys, "argv", ["gatk-sv-ploidy", "--help"])

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 0
    assert "Subcommands:" in capsys.readouterr().out


def test_main_rejects_unknown_subcommand(monkeypatch, capsys) -> None:
    monkeypatch.setattr(cli.sys, "argv", ["gatk-sv-ploidy", "bogus"])

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 1
    assert "unknown subcommand" in capsys.readouterr().out


def test_main_dispatches_to_subcommand(monkeypatch) -> None:
    called = {}

    def fake_main() -> None:
        called["main_called"] = True

    def fake_import(name: str):
        called["module_name"] = name
        return types.SimpleNamespace(main=fake_main)

    monkeypatch.setattr(cli.sys, "argv", ["gatk-sv-ploidy", "infer", "--input", "x"])
    monkeypatch.setattr("importlib.import_module", fake_import)

    cli.main()

    assert called == {
        "module_name": "gatk_sv_ploidy.infer",
        "main_called": True,
    }
    assert cli.sys.argv == ["gatk-sv-ploidy infer", "--input", "x"]
