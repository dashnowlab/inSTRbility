#!/usr/bin/env python
"""
Simple python wrapper for SSW library
Please put the path of libssw.so into LD_LIBRARY_PATH or pass it explicitly as a parameter
By Yongan Zhao (March 2016)
Revised by Mengyao Zhao on 2022-May-23
"""

import sys, os
import os.path as op
import argparse as ap
import ctypes as ct
import timeit as ti
import gzip
import math
from . import ssw_lib


def suppress_stderr():
    """Redirects stderr to /dev/null to suppress C++ error messages."""
    sys.stderr.flush()  # Flush any buffered stderr data
    stderr_fileno = sys.stderr.fileno()  # Get stderr file descriptor
    devnull = os.open(os.devnull, os.O_WRONLY)  # Open /dev/null for writing
    os.dup2(devnull, stderr_fileno)  # Redirect stderr to /dev/null
    os.close(devnull)  # Close descriptor

def restore_stderr(original_stderr):
    """Restores the original stderr after suppression."""
    sys.stderr.flush()
    os.dup2(original_stderr, sys.stderr.fileno())
    os.close(original_stderr)
    # sys.stderr = sys.__stderr__  # Restore original stderr


def read(sFile):
    """
    read a sequence file
    @param  sFile   sequence file
    """
    def read_one_fasta(f):
        """
        read a fasta file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        for l in f:
            if l.startswith('>'):
                if sSeq:
                    yield sId, sSeq, ''
                sId = l.strip()[1:].split()[0]
                sSeq = ''
            else:
                sSeq += l.strip()

        yield sId, sSeq, ''

    def read_one_fastq(f):
        """
        read a fastq file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        s3 = ''
        sQual = ''
        for l in f:
            sId = l.strip()[1:].split()[0]
            sSeq = f.readline().strip()
            s3 = f.readline().strip()
            sQual = f.readline().strip()

            yield sId, sSeq, sQual

# test if fasta or fastq
    bFasta = True
    ext = op.splitext(sFile)[1][1:].strip().lower()
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            l = f.readline()
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                sys.stderr.write('file format cannot be recognized\n')
                sys.exit()
    else:
        with open(sFile, 'r') as f:
            l = f.readline()
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                sys.stderr.write('file format cannot be recognized\n')
                sys.exit()

# read
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            if bFasta == True:
                for sId,sSeq,sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId,sSeq,sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual
    else:
        with open(sFile, 'r') as f:
            if bFasta == True:
                for sId,sSeq,sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId,sSeq,sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual


def to_int(seq, lEle, dEle2Int):
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i,ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n

    return num


def align_one(ssw, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    """
    align one pair of sequences
    @param  qProfile   query profile
    @param  rNum   number array for reference
    @param  nRLen   length of reference sequence
    @param  nFlag   alignment flag
    @param  nMaskLen   mask length
    """
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, nMaskLen)

    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)

    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)


def buildPath(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = 'MIDNSHP=X'
    sCigar = ''
    sQ = ''
    sA = ''
    sR = ''
    nQOff = nQryBeg
    nROff = nRefBeg
    for x in lCigar:
        n = x >> 4
        m = x & 15
        if m > 8:
            c = 'M'
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == 'M':
            sQ += q[nQOff : nQOff+n]
            sA += ''.join(['|' if q[nQOff+j] == r[nROff+j] else '*' for j in range(n)])
            sR += r[nROff : nROff+n]
            nQOff += n
            nROff += n
        elif c == 'I':
            sQ += q[nQOff : nQOff+n]
            sA += ' ' * n
            sR += '-' * n
            nQOff += n
        elif c == 'D':
            sQ += '-' * n
            sA += ' ' * n
            sR += r[nROff : nROff+n]
            nROff += n
    return sCigar, sQ, sA, sR


def include_substitution(sCigar, target, query, tbegin, qbegin):
    """
    Check if the CIGAR string includes substitutions
    @param  sCigar   CIGAR string
    @param  target   target sequence
    @param  query    query sequence
    """

    new_cigar = ""
    tpos = tbegin - 1  # convert to 0-based index
    qpos = qbegin - 1  # convert to 0-based index
    clen = ''; ctype = ''
    for c in sCigar:
        if c.isdigit():
            clen += c
        else:
            ctype = c
            if ctype == 'M' or ctype == '=':
                mlen = 0; slen = 0
                for i in range(int(clen)):
                    if target[tpos + i] == query[qpos + i]:
                        if slen > 0: new_cigar += str(slen) + 'X'
                        mlen += 1
                        slen = 0
                    else:
                        if mlen > 0: new_cigar += str(mlen) + 'M'
                        mlen = 0
                        slen += 1
                if mlen > 0: new_cigar += str(mlen) + 'M'
                if slen > 0: new_cigar += str(slen) + 'X'
                tpos += int(clen); qpos += int(clen)
            elif ctype == 'I':
                new_cigar += (str(clen) + 'I'); qpos += int(clen)
            elif ctype == 'D':
                new_cigar += (str(clen) + 'D'); tpos += int(clen)
            elif ctype == 'S': new_cigar += (str(clen) + 'S')
            clen = ''
    return new_cigar

def align_pair(target, query, nMatch=2, nMismatch=2, nOpen=3, nExt=1):
    lEle = []
    dRc = {} 
    dEle2Int = {}
    dInt2Ele = {}

    lEle = ['A', 'C', 'G', 'T', 'N']
    dRc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'T', 'c':'G', 'g':'C', 't':'A', 'N': 'N'} 
    for i,ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    nEleNum = len(lEle)
    lScore = [0 for i in range(nEleNum**2)]
    for i in range(nEleNum-1):
        for j in range(nEleNum-1):
            if lEle[i] == lEle[j]:
                lScore[i*nEleNum+j] = nMatch
            else:
                lScore[i*nEleNum+j] = -nMismatch


    # translate score matrix to ctypes
    mat = (len(lScore) * ct.c_int8) ()
    mat[:] = lScore
    # set flag
    nFlag = 2

    ssw = ssw_lib.CSsw(op.dirname(__file__))
    # iterate query sequence
    # build query profile
    # build rc query profile
    qNum = to_int(query, lEle, dEle2Int)
    qProfile = ssw.ssw_init(qNum, ct.c_int32(len(query)), mat, len(lEle), 2)
# build rc query profile
    query_rc = ''.join([dRc[x] for x in query[::-1]])
    qRcNum = to_int(query_rc, lEle, dEle2Int)
    qRcProfile = ssw.ssw_init(qRcNum, ct.c_int32(len(query_rc)), mat, len(lEle), 2)
# set mask len
    nMaskLen = len(query) // 2

    # iter target sequence
    rNum = to_int(target, lEle, dEle2Int)
    # format of res: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
    if nMaskLen < 15:
        original_stderr = os.dup(sys.stderr.fileno())
        suppress_stderr()
        res = align_one(ssw, qProfile, rNum, len(target), nOpen, nExt, nFlag, nMaskLen)
        restore_stderr(original_stderr)
    else:
        res = align_one(ssw, qProfile, rNum, len(target), nOpen, nExt, nFlag, nMaskLen)

# align rc query
    resRc = None
    if nMaskLen < 15:
        original_stderr = os.dup(sys.stderr.fileno())
        suppress_stderr()
        resRc = align_one(ssw, qRcProfile, rNum, len(target), nOpen, nExt, nFlag, nMaskLen)
        restore_stderr(original_stderr)
    else:
        resRc = align_one(ssw, qRcProfile, rNum, len(target), nOpen, nExt, nFlag, nMaskLen)

    # build cigar and trace back path
    strand = 0
    if resRc == None or res[0] >= resRc[0]:
        resPrint = res
        strand = 0
        sCigar, sQ, sA, sR = buildPath(query, target, res[4], res[2], res[8])
    else:
        resPrint = resRc
        strand = 1
        sCigar, sQ, sA, sR = buildPath(query_rc, target, resRc[4], resRc[2], resRc[8])
    
    alignment_score = resPrint[0]
    strand = '+' if strand == 0 else '-'
    target_begin = resPrint[2] + 1
    target_end = resPrint[3] + 1
    query_begin = resPrint[4] + 1
    query_end = resPrint[5] + 1

    # print('Alignment Score: {}'.format(alignment_score))
    # print('Strand: {}'.format(strand))
    # print('Target Begin: {}'.format(target_begin))
    # print('Target End: {}'.format(target_end))
    # print('Query Begin: {}'.format(query_begin))
    # print('Query End: {}'.format(query_end))
    # print('CIGAR: {}'.format(sCigar))

    ssw.init_destroy(qRcProfile)
    sCigar = include_substitution(sCigar, target, query, target_begin, query_begin)
    return (alignment_score, strand, target_begin, target_end, query_begin, query_end, sCigar)


if __name__ == '__main__':
    pass
