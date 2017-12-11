xiSPEC_ms_parser

mkdir ../uploads

virtualenv python_env

source python_env/bin/activate

pip install -r requirements.txt //pyteomics module has some changes so for the time being you need to checkout HEAD again afterwards

sqlite3 dbs/xiSPEC.db

CREATE TABLE "access_log" ( `id` INTEGER PRIMARY KEY AUTOINCREMENT, `ip` TEXT, `date` TEXT, `db_id` INTEGER );

CREATE TABLE "databases" ( `id` INTEGER PRIMARY KEY AUTOINCREMENT, `name` TEXT UNIQUE, `pass` TEXT, `share` TEXT UNIQUE, `ip` TEXT, `date` TEXT );
