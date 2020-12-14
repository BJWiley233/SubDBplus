#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 06:34:43 2020

@author: coyote
"""

import mysql.connector
from mysql.connector import errorcode


def connect(user, password, host):
    cnx = mysql.connector.connect(user=user, 
                                  password=password,
                                  host=host)
    return cnx
  
  
def connect_with_db(user, password, host, db_name):
    cnx = mysql.connector.connect(user=user, 
                                  password=password,
                                  host=host,
                                  database=db_name) 
    return cnx
    

def create_database(cursor, db_name):
    try:
        cursor.execute(
            "CREATE DATABASE IF NOT EXISTS {}".format(db_name))
    except mysql.connector.Error as err:
        print("Failed creating database: {}".format(err))

# private can only use here
def _drop_database(cursor, db_name):
    try:
        cursor.execute(
            "DROP DATABASE IF EXISTS {}".format(db_name))
    except mysql.connector.Error as err:
        print("Failed dropping database: {}".format(err))


def use_database(cnx, cursor, db_name):
    try:
        cursor.execute(
            "USE {}".format(db_name))
    except mysql.connector.Error as err:
        print("Database {} does not exist".format(db_name))
        if err.errno == errorcode.ER_BAD_DB_ERROR:
            create_database(cursor, db_name)
            print("Database {} created successfully".format(db_name))
            cnx.database = db_name
   
        
def drop_table(cursor, db_name, table):
    try:
        cursor.execute(
            "DROP TABLE IF EXISTS {}".format(table))
    except mysql.connector.Error as err:
        print("Something went wrong: {}".format(err))

