/* @source powerwater application
**
** 
** Application for large number of pairwise alignments without
** excessive filesystem traffic.
** Adapted from powerneedle by
** @author Copyright (C) Mayo Röttger (maroe001@uni-duesseldorf.de)
** 
**
**
** True Smith-Waterman best local alignment
** @author Copyright (C) Alan Bleasby (ableasby@hgmp.mrc.ac.uk)
** @@
**
** This program is free software; you can redistribute it and/or
** modify it under the terms of the GNU General Public License
** as published by the Free Software Foundation; either version 2
** of the License, or (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
******************************************************************************/

#include "emboss.h"




/*-
 * Copyright (c) 1991, 1993
 *	The Regents of the University of California.  All rights reserved.
 *
 * The folowing HEAPSORT code is derived from software contributed to
 * Berkeley by Ronnie Kon at Mindcraft Inc., Kevin Lew and Elmer Yglesias.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *	This product includes software developed by the University of
 *	California, Berkeley and its contributors.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * $DragonFly: src/lib/libc/stdlib/heapsort.c,v 1.5 2005/11/20 12:37:48 swildner Exp $
 *
 * @(#)heapsort.c	8.1 (Berkeley) 6/4/93
 */

#include <errno.h>
#include <stddef.h>
#include <stdlib.h>

/*
 * Swap two areas of size number of bytes.  Although qsort(3) permits random
 * blocks of memory to be sorted, sorting pointers is almost certainly the
 * common case (and, were it not, could easily be made so).  Regardless, it
 * isn't worth optimizing; the SWAP's get sped up by the cache, and pointer
 * arithmetic gets lost in the time required for comparison function calls.
 */
#define	SWAP(a, b, count, size, tmp) { \
    count = size; \
    do { \
        tmp = *a; \
        *a++ = *b; \
        *b++ = tmp; \
    } while (--count); \
}

/* Copy one block of size size to another. */
#define COPY(a, b, count, size, tmp1, tmp2) { \
    count = size; \
    tmp1 = a; \
    tmp2 = b; \
    do { \
        *tmp1++ = *tmp2++; \
    } while (--count); \
}

/*
 * Build the list into a heap, where a heap is defined such that for
 * the records K1 ... KN, Kj/2 >= Kj for 1 <= j/2 <= j <= N.
 *
 * There two cases.  If j == nmemb, select largest of Ki and Kj.  If
 * j < nmemb, select largest of Ki, Kj and Kj+1.
 */
#define CREATE(initval, nmemb, par_i, child_i, par, child, size, count, tmp) { \
    for (par_i = initval; (child_i = par_i * 2) <= nmemb; \
    par_i = child_i) { \
        child = base + child_i * size; \
        if (child_i < nmemb && compar(child, child + size) < 0) { \
            child += size; \
            ++child_i; \
        } \
        par = base + par_i * size; \
        if (compar(child, par) <= 0) \
            break; \
        SWAP(par, child, count, size, tmp); \
    } \
}

/*
 * Select the top of the heap and 'heapify'.  Since by far the most expensive
 * action is the call to the compar function, a considerable optimization
 * in the average case can be achieved due to the fact that k, the displaced
 * elememt, is ususally quite small, so it would be preferable to first
 * heapify, always maintaining the invariant that the larger child is copied
 * over its parent's record.
 *
 * Then, starting from the *bottom* of the heap, finding k's correct place,
 * again maintianing the invariant.  As a result of the invariant no element
 * is 'lost' when k is assigned its correct place in the heap.
 *
 * The time savings from this optimization are on the order of 15-20% for the
 * average case. See Knuth, Vol. 3, page 158, problem 18.
 *
 * XXX Don't break the #define SELECT line, below.  Reiser cpp gets upset.
 */
#define SELECT(par_i, child_i, nmemb, par, child, size, k, count, tmp1, tmp2) { \
    for (par_i = 1; (child_i = par_i * 2) <= nmemb; par_i = child_i) { \
        child = base + child_i * size; \
        if (child_i < nmemb && compar(child, child + size) < 0) { \
            child += size; \
            ++child_i; \
        } \
        par = base + par_i * size; \
        COPY(par, child, count, size, tmp1, tmp2); \
    } \
    for (;;) { \
        child_i = par_i; \
        par_i = child_i / 2; \
        child = base + child_i * size; \
        par = base + par_i * size; \
        if (child_i == 1 || compar(k, par) < 0) { \
            COPY(child, k, count, size, tmp1, tmp2); \
            break; \
        } \
        COPY(child, par, count, size, tmp1, tmp2); \
    } \
}

/*
 * Heapsort -- Knuth, Vol. 3, page 145.  Runs in O (N lg N), both average
 * and worst.  While heapsort is faster than the worst case of quicksort,
 * the BSD quicksort does median selection so that the chance of finding
 * a data set that will trigger the worst case is nonexistent.  Heapsort's
 * only advantage over quicksort is that it requires little additional memory.
 */
int
heapsort(void *vbase, size_t nmemb, size_t size,
    	 int (*compar)(const void *, const void *))
{
    int cnt, i, j, l;
    char tmp, *tmp1, *tmp2;
    char *base, *k, *p, *t;
    
    if (nmemb <= 1)
        return (0);
    
    if (!size) {
        errno = EINVAL;
        return (-1);
    }
    
    if ((k = malloc(size)) == NULL)
        return (-1);
    
    /*
     * Items are numbered from 1 to nmemb, so offset from size bytes
     * below the starting address.
     */
    base = (char *)vbase - size;
    
    for (l = nmemb / 2 + 1; --l;)
        CREATE(l, nmemb, i, j, t, p, size, cnt, tmp);
    
    /*
     * For each element of the heap, save the largest element into its
     * final slot, save the displaced element (k), then recreate the
     * heap.
     */
    while (nmemb > 1) {
        COPY(k, base + nmemb * size, cnt, size, tmp1, tmp2);
        COPY(base + nmemb * size, base + size, cnt, size, tmp1, tmp2);
        --nmemb;
        SELECT(i, j, nmemb, t, p, size, k, cnt, tmp1, tmp2);
    }
    free(k);
    return (0);
}






/* @func compare *************************************************************
**
** Compares two sequence names.
**
** @param [r] el1 [AjPSeq *] Pointer to sequence1
** @param [r] el2 [AjPSeq *] Pointer to sequence2
** @return [int] < 0, > 0 or 0 if name of sequence1 comes alphabetically before
** name of sequence2, after sequence2 or if the two sequences have identical
** names respectively.
** @@
******************************************************************************/

int compare(AjPSeq *el1, AjPSeq *el2)
{
    /* compare the two sequence names */
    return ajStrCmpS(ajSeqGetNameS(*el1),ajSeqGetNameS(*el2));
}

/* @func findRec *************************************************************
**
** Recursively finds a sequence by name in an alphabetically sorted sequence
** database and returns it.
** Because this is a pointer to the real internal sequence
** the caller must take care not to change the data in any way.
**
** @param [r] database [AjPSeqset] Sequence database
** @param [r] name [char *] Sequence name to search for
** @param [r] li [ajint] left boundary in database
** @param [r] re [ajint] right boundary in database
** @return [const AjPSeq] Pointer to the sequence in the database or NULL if
** sequence was not found
** @@
******************************************************************************/

const AjPSeq findRec(AjPSeqset database,char *name,ajint li,ajint re)
{    
    const AjPSeq seq;
    ajint mid;
    
    /* get sequence in the middle of boundaries */
    mid=(li+re)/2;
    seq=ajSeqsetGetseqSeq(database,mid);
    
    if (li<=re)
    {
        if (ajStrCmpC(ajSeqGetNameS(seq),name)>0) return findRec(database,name,li,mid-1); /* search left */
        else if (ajStrCmpC(ajSeqGetNameS(seq),name)<0) return findRec(database,name,mid+1,re);  /* search right */
        else return seq; /* sequence found */
    }
    else return NULL; /* sequence not found */
}

/* @func find ****************************************************************
 **
 ** Finds a sequence by name in a alphabetically sorted sequence database and
 ** returns it.
 ** Because this is a pointer to the real internal sequence
 ** the caller must take care not to change the data in any way.
 **
 ** @param [r] database [AjPSeqset] Sequence database
 ** @param [r] name [char *] Sequence name to search for
 ** @return [const AjPSeq] Pointer to the sequence in the database or NULL if
 ** sequence was not found
 ** @@
 *****************************************************************************/

const AjPSeq find(AjPSeqset database,char *name)
{
    return findRec(database,name,0,ajSeqsetGetSize(database)-1);
}






/* @prog powerwater ****************************************************************
**
** large number of Smith-Waterman local alignments
**
******************************************************************************/

int main(int argc, char **argv)
{
    AjPAlign align = NULL;
    AjPSeqset database =NULL;

    const AjPSeq finda = NULL;
    const AjPSeq findb = NULL;
    AjPSeq a = NULL;
    AjPSeq b = NULL;

    AjPStr m = NULL;
    AjPStr n = NULL;
    AjPStr ss = NULL;
    AjPStr fm = NULL;
    AjPStr fn = NULL;
    
    AjPFile identities = NULL;
    AjPFile pairs = NULL;
    
    AjBool error = ajFalse;
    AjBool dobrief = ajTrue;
    AjBool lastseq = ajFalse;

    const char  *p    = NULL;
    const char  *q    = NULL;
    char line[1024] = "";
    char seq1[512] = "";
    char seq2[512] = "";
    char c;
    
    ajint alignmentLength = 0;
    ajint start1 = 0;
    ajint start2 = 0;
    ajint end1 = 0;
    ajint end2 = 0;
    ajint k = 0;
    ajint *compass;
    
    ajuint mgaps  = 0;
    ajuint ngaps  = 0;
    ajuint olen;
    ajuint lena;
    ajuint lenb;
    ajuint i = 0;
    ajuint j = 1;
    ajuint linecount = 0;
    
    ajulong maxarr = 1000; 	/* arbitrary. realloc'd if needed */
    ajulong len;
    
    float *path;
    float **sub;
    float gapopen;
    float gapextend;
    float score;
    float id   = 0.;
    float sim  = 0.;
    float idx  = 0.;
    float simx = 0.;
    
    AjPMatrixf matrix;
    AjPSeqCvt cvt = 0;

    size_t stlen;



    /* Init everything */
    embInit("powerwater", argc, argv);
    
    ajFmtPrint("Database has been read.\n");
    
    matrix    = ajAcdGetMatrixf("datafile");
    database  = ajAcdGetSeqset("database");
    gapopen   = ajAcdGetFloat("gapopen");
    gapextend = ajAcdGetFloat("gapextend");
    dobrief   = ajAcdGetBoolean("brief");
    align     = ajAcdGetAlign("alignment");
    identities = ajAcdGetOutfile("identities");
    pairs     = ajAcdGetInfile("pairs");
    
    gapopen = ajRoundFloat(gapopen, 8);
    gapextend = ajRoundFloat(gapextend, 8);
    
    AJCNEW(path, maxarr);
    AJCNEW(compass, maxarr);
    
    m  = ajStrNew();
    n  = ajStrNew();
    ss = ajStrNew();
    
    sub = ajMatrixfGetMatrix(matrix);
    cvt = ajMatrixfGetCvt(matrix);
    
    /* Output Header for pairwise identity file */
    ajFmtPrintF(identities,"Sequence1 Sequence2 Identity Similarity Score Length #id #sim Start_Seq1 End_Seq1 Start_Seq2 End_Seq2\n");
    
    if (pairs!=NULL) /* if pairs file present */
    {
        /* Alphabetically sort database. This changes the sequence pointers of the database sequence set! */
        ajFmtPrint("Sorting database...\n");
        heapsort(database->Seq,ajSeqsetGetSize(database),sizeof(AjPSeq),(int(*)(const void*, const void*)) compare);
    }
    
    ajFmtPrint("Doing pairwise alignments...\n");    
    
    /*bis hier*/
    while (!lastseq) {
        error=ajFalse;
        
        /* Get next sequence pair */
        if (pairs!=NULL) /* Align pairs in pairsfile */
        {
            /* read next line from file */
            line[0]='\0';
            i=0;
            while ( (c = fgetc(ajFileGetFileptr(pairs)) ) != EOF && c != '\n') line[i++]=c;
            line[i]='\0';
            linecount++;
            if (feof(ajFileGetFileptr(pairs))) lastseq=ajTrue;
            
            if (strcmp(line,"")==0) error=ajTrue; /* not process empty line then */
            else
            {
                /* load first string */
                i=0;
                j=0;
                while ((c=line[i])!='\0' && line[i]!='\t' && line[i]!=' ')
                {
                    seq1[j++]=c;
                    i++;
                }
                seq1[j]='\0';
                
                /* load second string */
                i++;
                j=0;
                if (c!='\0') {
                    while ((c=line[i])!='\0' && line[i]!='\t' && line[i]!=' ')
                    {
                        seq2[j++]=c;
                        i++;
                    }
                }
                seq2[j]='\0';
                
                if (strcmp(seq1,"")==0)
                {
                    ajFmtPrint("Error! Could not read sequence 1 from pairs file line %d! Sequence 2 is '%s'. Skipping this pair.\n",linecount,seq2);
                    error=ajTrue;
                }
                
                if (strcmp(seq2,"")==0)
                {
                    ajFmtPrint("Error! Could not read sequence 2 from pairs file line %d! Sequence 1 is '%s'. Skipping this pair.\n",linecount,seq1);
                    error=ajTrue;
                }
            
                
                if (!error)
                {
                    /* find sequences in alphabetically sorted sequence set */
                    if ((finda=find(database,seq1))==NULL)
                    {
                        ajFmtPrint("Error! Could not find '%s' in database! Skipping pair '%s'-'%s' from pairs file line %d.\n",seq1,seq1,seq2,linecount);
                        error=ajTrue;
                    }
                    
                    if ((findb=find(database,seq2))==NULL)
                    {
                        ajFmtPrint("Error! Could not find '%s' in database! Skipping pair '%s'-'%s' from pairs file line %d.\n",seq2,seq1,seq2,linecount);
                        error=ajTrue;
                    }
                    
                    if (!error)
                    {
                        a=ajSeqNewSeq(finda);
                        b=ajSeqNewSeq(findb);
                    }
                }
            }
        }
        else /* Align all sequence pairs */
        {
            /* take current sequence pair i,j */
            a=ajSeqNewSeq(ajSeqsetGetseqSeq(database,i));
            b=ajSeqNewSeq(ajSeqsetGetseqSeq(database,j++));
            
            /* update indices */
            if (j==ajSeqsetGetSize(database))
            {
                if (i==j-2) lastseq=ajTrue;
                else j=(++i)+1;
            }            
        }
        
        
	     /* Do Smith Waterman pairwise local alignment */
	    if(!error)
	    {
			ajSeqTrim(a);
            ajSeqTrim(b);
            
            lena = ajSeqGetLen(a);
            lenb = ajSeqGetLen(b);
            
            mgaps=0;
            ngaps=0;
            
			if(lenb > (ULONG_MAX/(ajulong)(lena+1)))
			{
				//ajDie("Sequences too big. Try 'matcher' or 'supermatcher'");
				ajFmtPrint("Error in aligning of '%s'-'%s', pairs file line %d: Sequences too big. Try 'stretcher' or 'supermatcher' for this sequence pair!\n",ajSeqGetNameC(a),ajSeqGetNameC(b),linecount);
                error=ajTrue;
            }
            else
            {
				len = lena*lenb;
				
                if(len>maxarr)
                {
                    stlen = (size_t) len;
                    AJCRESIZETRY(path,stlen);
                    if(!path)
                    {
                        /* ajDie("Sequences too big. Try 'stretcher'"); */
                        ajFmtPrint("Error in aligning '%s'-'%s', pairs file line %d: Sequences too big. Try 'stretcher' for this sequence pair!\n",ajSeqGetNameC(a),ajSeqGetNameC(b),linecount);
                        error=ajTrue;
                    }
                    AJCRESIZETRY(compass,stlen);
                    if(!compass)
                    {
                        /* ajDie("Sequences too big. Try 'stretcher'"); */
                        ajFmtPrint("Error in aligning '%s'-'%s': Sequences too big. Try 'stretcher' for this sequence pair!\n",ajSeqGetNameC(a),ajSeqGetNameC(b),linecount);
                        error=ajTrue;
                    }
                    maxarr=len;
                }
            }

		//beginb=ajSeqGetBegin(b)+ajSeqGetOffset(b);
			if (!error)
            {	
				p = ajSeqGetSeqC(a);
				q = ajSeqGetSeqC(b);
	
				ajStrAssignC(&m,"");
				ajStrAssignC(&n,"");
	
				score = embAlignPathCalcSW(p,q,lena,lenb,gapopen,gapextend,path,sub,cvt,
				   compass,ajFalse);//show > ajFalse
	
				/*score=embAlignScoreSWMatrix(path,compass,gapopen,gapextend,a,b,lena,
					lenb,sub,cvt,&start1,&start2);*/
	
				embAlignWalkSWMatrix(path,compass,gapopen,gapextend,a,b,&m,&n,
				     lena,lenb,&start1,&start2);
	
		//~ ajDebug("ReportLocal call start1:%d begina:%d start2:%d beginb:%d\n",
			//~ start1, begina, start2, beginb);
		//~ 
				if(!dobrief)
				{
					embAlignReportLocal(align, a, b, m, n,
					    start1, start2,
					    gapopen, gapextend,
					    score, matrix,  ajSeqGetOffset(a), ajSeqGetOffset(b));
					    
					ajAlignWrite(align);
					ajAlignReset(align);
					
					/*    
					embAlignCalcSimilarity(m,n,sub,cvt,lena,lenb,&id,&sim,&idx,
							   &simx);
					ajFmtPrintS(&tmpstr,"Longest_Identity = %5.2f%%\n",
						   id);
					ajFmtPrintAppS(&tmpstr,"Longest_Similarity = %5.2f%%\n",
						   sim);
					ajFmtPrintAppS(&tmpstr,"Shortest_Identity = %5.2f%%\n",
						   idx);
					ajFmtPrintAppS(&tmpstr,"Shortest_Similarity = %5.2f%%",
						   simx);
					ajAlignSetSubHeaderApp(align, tmpstr);
					*/
				}
				
				/* Calculate identity and similarity */
                embAlignCalcSimilarity(m,n,sub,cvt,lena,lenb,&id,&sim,&idx,
                                       &simx);
                
                
                /* Count gaps in alignment string */
                ajStrAssignS(&fm,m);
                ajStrAssignS(&fn,n);
                ajStrFmtUpper(&fm);
                ajStrFmtUpper(&fn);
                
                p = ajStrGetPtr(fm);
                q = ajStrGetPtr(fn);
                olen = (ajint) strlen(p);
                
                for(k=0;k<olen;++k)
                {
                    if(p[k] =='.') ++mgaps;
                    if(q[k] =='.') ++ngaps;
                }
                
                /* Calculate alignment length - wrong b/c local alignment -fix it */
                end1=start2+ajSeqGetLen(a)+mgaps;
                end2=start1+ajSeqGetLen(b)+ngaps;
                alignmentLength=end1>end2?end1:end2;
                
                
                
                /* Adjust percent identity and similarity to the basis of alignment length */
                id=idx*(ajSeqGetLen(a)>ajSeqGetLen(b)?ajSeqGetLen(a):ajSeqGetLen(b))/alignmentLength;
                sim=simx*(ajSeqGetLen(a)>ajSeqGetLen(b)?ajSeqGetLen(a):ajSeqGetLen(b))/alignmentLength;
                
                
                
                /* Output pairwise identities */
                ajFmtPrintF(identities,"%s %s %.1f %.1f %.1f %d %.0f %.0f %d %d %d %d\n",ajSeqGetNameC(a),ajSeqGetNameC(b),id,sim,score,alignmentLength,id*alignmentLength/100,sim*alignmentLength/100,start1,end1,start2,end2);
            }
            else
            {
                /* Output pairwise '-' for this sequence pair */
                ajFmtPrintF(identities,"%s %s - - - - - - - - - -\n",ajSeqGetNameC(a),ajSeqGetNameC(b));
            }
        }
	
	    
    /* ab hier */
        /* destroy sequence copies */
        ajSeqDel(&a);
        ajSeqDel(&b);
	}
    /* Destroy everything */
    ajSeqsetDel(&database);
    
    ajAlignClose(align);
    ajAlignDel(&align);
    
    ajFileClose(&identities);
    if (pairs!=NULL) ajFileClose(&pairs);
    
    AJFREE(compass);
    AJFREE(path);
    
    ajStrDel(&n);
    ajStrDel(&m);
    ajStrDel(&ss);
    ajStrDel(&fm);
    ajStrDel(&fn);
    
    ajFmtPrint("Done.\n");
    
    embExit();
    
    return 0;
}
