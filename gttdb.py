"""Backward-compatible wrapper for the renamed `gtdbkit` CLI."""

from gtdbkit import main


if __name__ == "__main__":
    raise SystemExit(main())
