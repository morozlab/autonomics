#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      mat
#
# Created:     25/07/2012
# Copyright:   (c) mat 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import imaplib
import sys
import os

def main():
    c = imaplib.IMAP4_SSL('imap.gmail.com', 993)
    c.login('morozgenomics', 'Whitney2012')
    c.list()
    typ, data = c.select("Inbox")
    print("Messages in mailbox: "  + str(data[0]))
    typ, msg_ids = c.search(None, '(SUBJECT "PBS JOB 23504899.torx.ufhpc" BODY "exceeded MEM usage hard limit")')
    print(msg_ids)
    print(str(len(msg_ids)))
    if(msg_ids[0] == ''):
        print('moo')




if __name__ == '__main__':
    main()
