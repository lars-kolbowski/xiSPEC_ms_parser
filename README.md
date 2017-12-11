xiSPEC_ms_parser

mkdir dbs

mkdir dbs/tmp

mkdir dbs/saved

mkdir ../uploads

virtualenv python_env

source python_env/bin/activate

pip install -r requirements.txt


sqlite3 dbs/xiSPEC.db

CREATE TABLE "access_log" ( `id` INTEGER PRIMARY KEY AUTOINCREMENT, `ip` TEXT, `date` TEXT, `db_id` INTEGER );

CREATE TABLE "databases" ( `id` INTEGER PRIMARY KEY AUTOINCREMENT, `name` TEXT UNIQUE, `pass` TEXT, `share` TEXT UNIQUE, `ip` TEXT, `date` TEXT );
