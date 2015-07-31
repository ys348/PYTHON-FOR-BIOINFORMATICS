
#!/usr/bin/python

# Import modules for CGI handling 

import cgi, cgitb 
import sqlite3

print "Content-type:text/html\r\n\r\n"

print "<html>"

print "<head>"

print "<title>Gene Query</title>"

print "</head>"

print "<body>"

form = cgi.FieldStorage() 

# Get data from fields
option = form['option'].value
query = form['query'].value

# create a new database called 'final.db'
connection=sqlite3.connect('final.db') 
cursor=connection.cursor()  

if option == "gene":
    # generate the picture to display 
    cursor.execute("SELECT * FROM CuffDiff WHERE test_id='%s'" % query) #query
    result = cursor.fetchall()
    print result
else:
    #query
    cursor.execute('select ...') 
    result = cursor.fetchall()

print "<form action=\"/cgi-bin/handle_query.py\">"

print "<label>option:</label>"

print "<select name=\"option\">"
print "<option value=\"gene\">GENE</option>"
print "<option value=\"keyword\">Keyword</option>"
print "</select>"

print "<input type=\"text\" name=\"query\" value=\"\">"
print "<input type=\"submit\" value=\"Submit\" />"
print "</form>"

if option == "gene":
    print "Draw a picture here"
else:
    print "keyword"

print "</body>"

print "</html>"