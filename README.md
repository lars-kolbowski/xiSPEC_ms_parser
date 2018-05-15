# xiSPEC_ms_parser

Back-end parser for xiSPEC mass spectrometry visualization tool.


### Requirements:
python2.7

virtualenv

sqlite3

### Installation

Clone git repository into your web-server directory (e.g. /var/www/html):

```git clone https://github.com/Rappsilber-Laboratory/xiSPEC_ms_parser.git```

cd into the repository:

```cd xiSPEC_ms_parser```

Create uploads directory (permissions?):

```mkdir ../uploads```

Change owner of uploads directory to www-data:

```sudo chown www-data:www-data ../uploads/```

Change owner of log directory to www-data:

```sudo chown www-data:www-data log```

Change owner of dbs directory (and sub directories) to www-data:

```sudo chown -R www-data:www-data dbs```



Create python virtualenv:

```virtualenv --no-site-packages python_env```

Activate virtualenv:

```source python_env/bin/activate```

Install dependencies:

```pip install -r requirements.txt```

pyteomics module has some changes so for the time being you need to reset to HEAD again after installing the dependencies

```git reset --hard HEAD```
