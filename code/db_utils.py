#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    Utilities for database manipulation (currently using SQLite database)
    

Explanation:

    

Example:

    

'''

import os
import sqlite3


# connect to SQLite database
def connect_db(db_file, host_name='', user_name='', db_name='', unix_socket_loc=''):
    db_file = os.path.abspath(db_file)
#     # TODO use relativeToAbsolutePath helper function above
    print 'using database file: ' + db_file
    if(not os.path.exists(os.path.dirname(db_file))):
        os.makedirs(os.path.dirname(db_file))
    return sqlite3.connect(db_file)

# gets a list of all tables in the given database
def get_table_names(db_cursor):
    c = db_cursor
    # get existing list of tables
    c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
    tables= c.fetchall()
    table_names = [t[0] for t in tables]
    return table_names

def getColumnNamesOfTable(db_connection_cursor, tablename):
    columnsInfo = getColumnsOfTable(db_connection_cursor, tablename)
    column_names = [c[1] for c in columnsInfo]
    return column_names

#returns an array with info about each column (column[1] = name, column[2] = type)
def getColumnsOfTable(db_connection_cursor, tablename):
    try:
        db_connection_cursor.execute('pragma table_info({table});'.format(table=tablename))
        columns = db_connection_cursor.fetchall();
    except sqlite3.OperationalError:
        print 'error on tablename: ' + str(tablename)
        raise
    return columns

def is_region_dataset(dataset_name):
    return dataset_name.endswith('_r')
