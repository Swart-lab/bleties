#!/usr/bin/env python3

"""Variables shared between all Bleties scripts"""

class SharedValues():
    # Lists of reference- and query-consuming operations
    REFCONSUMING = ['M', 'D', 'N', '=', 'X']
    QUERYCONSUMING = ['M', 'I', 'S', '=', 'X']
    GFF3COLUMNS = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
