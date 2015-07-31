import sqlite3

def create_database():
	# create a new database called 'final.db'
    connection=sqlite3.connect('final.db') 
    cursor=connection.cursor()

	# open the file
    f=open('gene_exp.diff','r') 
    d=[i.strip().split('\t') for i in f.readlines()]
    f.close()
	# delete the header
    d=d[1:] 

    # create cuffDiff table
    cursor.execute('create table CuffDiff(test_id text, gene text, locus text, \
        sample_1 text, sample_2 text, status text, value_1 float, value_2 float, \
        ln_fc float, test_stat float, p_value float, significant text)') 

    #insert all data sets from 'gene_exp' into this new table
    cursor.executemany('INSERT INTO CuffDiff VALUES(?,?,?,?,?,?,?,?,?,?,?,?)',d)
    connection.commit()
    connection.close()
    #f = open('GeneTable', 'r')
    #TODO:truncate transcript (remove .X) to get gene name
    
    #TODO:create gene table
    #cursor.execute('create table ...')
    
    #TODO:insert records into database
    #cursor.executemany('insert ...')
    #cursor.commit()
    
if __name__ == "__main__":
    create_database() 