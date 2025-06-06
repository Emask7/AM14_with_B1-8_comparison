eh = ExperimentHub()
query(eh , 'msigdb')

msigdb.mm = getMsigdb(org = 'mm', id = 'SYM')
msigdb.mm
listCollections(msigdb.mm)
listSubCollections(msigdb.mm)
