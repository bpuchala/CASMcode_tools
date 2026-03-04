import os
import pathlib
import sqlite3
import time
import typing

root = pathlib.Path(os.environ["HOME"]) / ".config" / "casm"


def get_casm_config_dir():
    """Get the ~/.config/casm directory."""
    return root


def get_user_cache_path(cache: str) -> pathlib.Path:
    """Get the path to the sqlite database file for a given database name."""
    db_dir = get_casm_config_dir() / "sqlite_cache"
    db_dir.mkdir(parents=True, exist_ok=True)
    return db_dir / f"{cache}.db"


def get_local_cache_path(
    dir: typing.Union[pathlib.Path, str], cache: str
) -> pathlib.Path:
    """Get the path to the directory for a given cache name."""
    db_dir = dir / "sqlite_cache"
    db_dir.mkdir(parents=True, exist_ok=True)
    return db_dir / f"{cache}.db"


def init_db(db_path: str | pathlib.Path):
    """Create the sqlite database and table if it doesn't exist."""
    db_path = pathlib.Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS cache (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL,
            created_at REAL
        )
        """)
    conn.commit()
    conn.close()


def store_kv(db_path: str | pathlib.Path, key: str, value: str):
    """Store or replace a key/value pair into the DB."""
    init_db(db_path)
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute(
        "INSERT OR REPLACE INTO cache (key, value, created_at) VALUES (?, ?, ?)",
        (key, value, time.time()),
    )
    conn.commit()
    conn.close()


def get_value(db_path: str | pathlib.Path, key: str) -> str | None:
    """Retrieve the value for a given key, or None if not found."""
    init_db(db_path)
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()
    cur.execute("SELECT value FROM cache WHERE key = ?", (key,))
    row = cur.fetchone()
    conn.close()
    return row[0] if row else None


class UserCache:
    """Use simple key/value caches stored in a sqlite database in the user's home
    directory.

    Key/value pairs stored for cache=="foo" will be stored in a sqlite database file
    located at `~/.config/casm/sqlite_cache/foo.db`.
    """

    def __init__(self):
        pass

    def store(self, cache: str, key: str, value: str):
        """Store a key/value pair for a given request index."""
        db_path = get_user_cache_path(cache)
        store_kv(db_path, key, value)

    def get(self, cache: str, key: str) -> str | None:
        """Get the value for a given key and request index, or None if not found."""
        db_path = get_user_cache_path(cache)
        return get_value(db_path, key)


class LocalCache:
    """A simple cache that stores key/value pairs in sqlite databases.

    Key/value pairs stored for cache=="foo" will be stored in a sqlite database file
    located at `{dir}/sqlite_cache/foo.db`, where `dir` is the directory specified when
    initializing the LocalCache.
    """

    def __init__(self, dir: typing.Union[pathlib.Path, str]):
        self.dir = pathlib.Path(dir)
        self.dir.mkdir(parents=True, exist_ok=True)

    def store(self, cache: str, key: str, value: str):
        """Store a key/value pair for a given request index."""
        db_path = get_local_cache_path(self.dir, cache)
        store_kv(db_path, key, value)

    def get(self, cache: str, key: str) -> str | None:
        """Get the value for a given key and request index, or None if not found."""
        db_path = get_local_cache_path(self.dir, cache)
        return get_value(db_path, key)
