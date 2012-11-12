/* Module for pairwise sequence alignment in the LingPy library.
 */

# include <Python.h>
# include <boost/python.hpp>
# include <boost/python/dict.hpp>
# include <boost/python/list.hpp>
# include <boost/python/object.hpp>
# include <string>
# include <iostream>
# include <vector>

using namespace boost::python;

float _global(
        int *,
        int,
        int *,
        int,
        float **,
        float,
        int *,
        int *
        );

float _local(
        int *,
        int,
        int *,
        int,
        float **,
        float,
        int *,
        int *
        );

float _overlap(
        int *,
        int,
        int *,
        int,
        float **,
        float,
        int *,
        int *
        );

float _repeats(
        int *,
        int,
        int *,
        int,
        float **,
        float,
        int *,
        int *
        );

float _dialign(
        int *,
        int,
        int *,
        int,
        float **,
        float,
        int *,
        int *
        );

tuple align_pairwise(
        list,
        list,
        list,
        list,
        list,
        list,
        std::string,
        std::string,
        dict,
        float,
        float,
        std::string
        );

list align_sequences_pairwise(
        list, 
        list, 
        list,
        list, 
        dict,
        float,
        float,
        std::string
        );

list align_sequence_pairs(
        list, 
        list, 
        list,
        list, 
        dict, 
        float,
        float,
        std::string
        );

dict random_align_sequence_pairs(
        list,
        list,
        list,
        list,
        dict,
        float,
        float,
        std::string,
        int
        );

float edit_dist(
        list seqA,
        list seqB
        );

/* The basic functions for pairwise sequence alignment */

/* Function _global aligns sequences globally, following the idea of site
 * restrictions, which cannot insert gaps in specific sites of the other
 * sequences. */
float _global(
        int * arrayA,
        int lA,
        int * arrayB,
        int lB,
        float ** scorer,
        float scale,
        int * almA,
        int * almB
        )
{    
    int i,j;
    float gapA,gapB,match;

    // initialize matrix and traceback
    float matrix[lB+1][lA+1];
    int traceback[lB+1][lA+1];

    matrix[0][0] = 0.0;
    traceback[0][0] = 1;

    for(i=1;i<lA+1;i++)
    {
        matrix[0][i] = matrix[0][i-1] + scorer[0][i] * scale;
        traceback[0][i] = 2;
    }
    for(i=1;i<lB+1;i++)
    {
        matrix[i][0] = matrix[i-1][0] + scorer[i][0] * scale;
        traceback[i][0] = 3;
    }

    for(i=1;i<lB+1;i++)
    {
        for(j=1;j<lA+1;j++)
        {
            // calculate the costs for gapping seqB
            if(arrayB[i-1] < 0 && arrayA[j-1] > 0 && j != lA)
            {
                gapA = matrix[i-1][j] - 1000000;
            }
            else if(traceback[i-1][j] == 3)
            {
                gapA = matrix[i-1][j] + scorer[i][0] * scale;
            }
            else
            {
                gapA = matrix[i-1][j] + scorer[i][0];
            }
            
            // calculate the costs for gapping seqA
            if(arrayA[j-1] < 0 && arrayB[i-1] > 0 && i != lB)
            {
                gapB = matrix[i][j-1] - 1000000;
            }
            else if(traceback[i][j-1] == 2)
            {
                gapB = matrix[i][j-1] + scorer[0][j] * scale;
            }
            else
            {
                gapB = matrix[i][j-1] + scorer[0][j];
            }

            // calculate the costs for a match
            match = matrix[i-1][j-1] + scorer[i][j];

            // compare the costs
            if(gapA > match && gapA >= gapB)
            {
                matrix[i][j] = gapA;
                traceback[i][j] = 3;
            }
            else if(match >= gapB)
            {
                matrix[i][j] = match;
                traceback[i][j] = 1;
            }
            else
            {
                matrix[i][j] = gapB;
                traceback[i][j] = 2;
            }
        }
    }
    // calculate the similarity score
    float sim = matrix[lB][lA];
    
    // carry out the traceback
    i--;
    j--;

    while(i > 0 || j > 0)
    {
        if(traceback[i][j] == 3)
        {
            almA[j] = almA[j] + 1;
            i--;
        }
        else if(traceback[i][j] == 1)
        {
            i--;
            j--;
        }
        else
        {
            almB[i] = almB[i] + 1;
            j--;
        }
    }
    return sim;
}

float _local(
        int * arrayA,
        int lA,
        int * arrayB,
        int lB,
        float ** scorer,
        float scale,
        int * almA,
        int * almB
        )
{    
    int i,j,k;
    float gapA,gapB,match;
    float null;
    //@check float max;
    int imax,jmax;
    float max_score = 0.0;


    // initialize matrix and traceback
    float matrix[lB+1][lA+1];
    int traceback[lB+1][lA+1];

    matrix[0][0] = 0.0;
    traceback[0][0] = 0;

    for(i=1;i<lA+1;i++)
    {
        matrix[0][i] = 0;
        traceback[0][i] = 0;
    }
    for(i=1;i<lB+1;i++)
    {
        matrix[i][0] = 0;
        traceback[i][0] = 0;
    }

    for(i=1;i<lB+1;i++)
    {
        for(j=1;j<lA+1;j++)
        {
            null = 0.0;

            // calculate the costs for gapping seqB
            if(arrayB[i-1] < 0 && arrayA[j-1] > 0 && j != lA)
            {
                gapA = matrix[i-1][j] - 1000000;
                null = -1000000;
            }
            else if(traceback[i-1][j] == 3)
            {
                gapA = matrix[i-1][j] + scorer[i][0] * scale;
            }
            else
            {
                gapA = matrix[i-1][j] + scorer[i][0];
            }
            
            // calculate the costs for gapping seqA
            if(arrayA[j-1] < 0 && arrayB[i-1] > 0 && i != lB)
            {
                gapB = matrix[i][j-1] - 1000000;
                null = -1000000;
            }
            else if(traceback[i][j-1] == 2)
            {
                gapB = matrix[i][j-1] + scorer[0][j] * scale;
            }
            else
            {
                gapB = matrix[i][j-1] + scorer[0][j];
            }

            // calculate the costs for a match
            match = matrix[i-1][j-1] + scorer[i][j];

            // compare the costs
            if(gapA >= match && gapA >= gapB && gapA >= null )
            {
                matrix[i][j] = gapA;
                traceback[i][j] = 3;
            }
            else if(match >= gapB && match >= null)
            {
                matrix[i][j] = match;
                traceback[i][j] = 1;
            }
            else if(gapB > null)
            {
                matrix[i][j] = gapB;
                traceback[i][j] = 2;
            }
            else
            {
                matrix[i][j] = null;
                traceback[i][j] = 0;
            }

            // check for maximum score
            if(matrix[i][j] >= max_score)
            {
                max_score = matrix[i][j];
                imax = i;
                jmax = j;
            }
        }
    }
    // calculate the similarity score
    float sim = matrix[imax][jmax];
    
    // get the maximum values
    i = imax;
    j = jmax;

    // mark the upper boundaries of the local alignment
    for(k=j;k<lA;k++)
    {
        almA[k] = -1;
    }
    for(k=i;k<lB;k++)
    {
        almB[k] = -1;
    }

    // carry out the traceback
    while(traceback[i][j] != 0)
    {
        if(traceback[i][j] == 3)
        {
            almA[j] = almA[j] + 1;
            i--;
        }
        else if(traceback[i][j] == 1)
        {
            i--;
            j--;
        }
        else if(traceback[i][j] == 2)
        {
            almB[i] = almB[i] + 1;
            j--;
        }
        else
        {
            break;
        }
    }

    // mark the lower boundaries
    for(k=0;k<j;k++)
    {
        almA[k] = -1;
    }
    for(k=0;k<i;k++)
    {
        almB[k] = -1;
    }
    
    return sim;
}

float _overlap(
        int * arrayA,
        int lA,
        int * arrayB,
        int lB,
        float ** scorer,
        float scale,
        int * almA,
        int * almB
        )
{    
    int i,j;
    float gapA,gapB,match;

    // initialize matrix and traceback
    float matrix[lB+1][lA+1];
    int traceback[lB+1][lA+1];

    matrix[0][0] = 0.0;
    traceback[0][0] = 1;

    for(i=1;i<lA+1;i++)
    {
        matrix[0][i] = 0.0;
        traceback[0][i] = 2;
    }
    for(i=1;i<lB+1;i++)
    {
        matrix[i][0] = 0.0;
        traceback[i][0] = 3;
    }

    for(i=1;i<lB+1;i++)
    {
        for(j=1;j<lA+1;j++)
        {
            // calculate the costs for gapping seqB
            if(arrayB[i-1] < 0 && arrayA[j-1] > 0 && j != lA)
            {
                gapA = matrix[i-1][j] - 1000000;
            }
            else if(j == lA)
            {
                gapA = matrix[i-1][j];
            }
            else if(traceback[i-1][j] == 3)
            {
                gapA = matrix[i-1][j] + scorer[i][0] * scale;
            }
            else
            {
                gapA = matrix[i-1][j] + scorer[i][0];
            }
            
            // calculate the costs for gapping seqA
            if(arrayA[j-1] < 0 && arrayB[i-1] > 0 && i != lB)
            {
                gapB = matrix[i][j-1] - 1000000;
            }
            else if(i == lB)
            {
                gapB = matrix[i][j-1];
            }
            else if(traceback[i][j-1] == 2)
            {
                gapB = matrix[i][j-1] + scorer[0][j] * scale;
            }
            else
            {
                gapB = matrix[i][j-1] + scorer[0][j];
            }

            // calculate the costs for a match
            match = matrix[i-1][j-1] + scorer[i][j];

            // compare the costs
            if(gapA > match && gapA >= gapB)
            {
                matrix[i][j] = gapA;
                traceback[i][j] = 3;
            }
            else if(match >= gapB)
            {
                matrix[i][j] = match;
                traceback[i][j] = 1;
            }
            else
            {
                matrix[i][j] = gapB;
                traceback[i][j] = 2;
            }
        }
    }
    // calculate the similarity score
    float sim = matrix[lB][lA];
    
    // carry out the traceback
    i--;
    j--;

    while(i > 0 || j > 0)
    {
        if(traceback[i][j] == 3)
        {
            almA[j] = almA[j] + 1;
            i--;
        }
        else if(traceback[i][j] == 1)
        {
            i--;
            j--;
        }
        else
        {
            almB[i] = almB[i] + 1;
            j--;
        }
    }
    return sim;
}

float _repeats(
        int * arrayA,
        int lA,
        int * arrayB,
        int lB,
        float ** scorer,
        float scale,
        int * almA,
        int * almB
        )
{    
    int i,j,k,l,o;
    float gapA,gapB,match;
    float null;
    //@check float max;
    //@check int imax,jmax;
    //@check float max_score = 0.0;
    bool breakit = false;
    int skipA,skipB,skip;


    // initialize matrix and traceback
    float matrix[lB+1][lA+1];
    int traceback[lB+1][lA+1];

    matrix[0][0] = 0.0;
    traceback[0][0] = 1;

    for(i=1;i<lA+1;i++)
    {
        matrix[0][i] = 0;
        traceback[0][i] = 2;
    }
    for(i=1;i<lB+1;i++)
    {
        matrix[i][0] = 0;
        traceback[i][0] = 3;
    }

    for(i=1;i<lB+1;i++)
    {
        for(j=1;j<lA+1;j++)
        {
            // calculate the costs for gapping seqB
            if(arrayB[i-1] < 0 && arrayA[j-1] > 0 && j != lA)
            {
                gapA = matrix[i-1][j] - 1000000;
                null = -1000000;
            }
            else if(traceback[i-1][j] == 3)
            {
                gapA = matrix[i-1][j] + scorer[i][0] * scale;
                null = 0.0;
            }
            else
            {
                gapA = matrix[i-1][j] + scorer[i][0];
                null = 0.0;
            }
            
            // calculate the costs for gapping seqA
            if(arrayA[j-1] < 0 && arrayB[i-1] > 0 && i != lB)
            {
                gapB = matrix[i][j-1] - 1000000;
                null = -1000000;
            }
            else if(traceback[i][j-1] == 2)
            {
                gapB = matrix[i][j-1] + scorer[0][j] * scale;
                null = 0.0;
            }
            else
            {
                gapB = matrix[i][j-1] + scorer[0][j];
                null = 0.0;
            }

            // calculate the costs for a match
            match = matrix[i-1][j-1] + scorer[i][j];

            // compare the costs
            if(gapA >= match && gapA >= gapB && gapA >= null )
            {
                matrix[i][j] = gapA;
                traceback[i][j] = 3;
            }
            else if(match > gapB && match > null)
            {
                matrix[i][j] = match;
                traceback[i][j] = 1;
            }
            else if(gapB > null)
            {
                matrix[i][j] = gapB;
                traceback[i][j] = 2;
            }
            else
            {
                matrix[i][j] = null;
                traceback[i][j] = 0;
            }
        }
    }
    // calculate the similarity score
    float sim = matrix[lB][lA];

    i--;
    j--;

    // carry out the traceback
    while(i > 0 || j > 0)
    {
        if(traceback[i][j] == 3)
        {
            almA[j] = almA[j] + 1;
            i--;
        }
        else if(traceback[i][j] == 1)
        {
            i--;
            j--;
        }
        else if(traceback[i][j] == 2)
        {
            almB[i] = almB[i] + 1;
            j--;
        }
        else
        {
            for(k=j-1;k>=0;k--)
            {
                for(l=i-1;l>=0;l--)
                {
                    if(traceback[l][k] == 1);
                    breakit = true;
                    break;
                }
                if(breakit == true)
                {
                    breakit = false;
                    break;
                }
            }
            skip = (j-k)+(i-l);
            skipA = j-k;
            skipB = i-l;
            for(o=0;o<skip-skipA;o++)
            {
                almA[j-skipA] = almA[j-skipA] + 1; //(j-skipA,0);
            }
            for(o=0;o<skip-skipB;o++)
            {
                almB[i] = almB[i] + 1; //.insert(i,0);
            }
            i = l;
            j = k;
            sim = sim + matrix[i][j];
        }
    }
    
    return sim;
}

float _dialign(
        int * arrayA,
        int lA,
        int * arrayB,
        int lB,
        float ** scorer,
        float scale,
        int * almA,
        int * almB
        )
{

    list subA,subB;

    int i,j,k,l;
    int minimum;
    
    float old_score;
    float new_score;
    int old_length;
    int new_length;

    float scoreA,scoreB;
    float max_score;
    float sim;

    /* create the matrix */
    float matrix[lB+1][lA+1];
    int traceback[lB+1][lA+1];
 
    /* initialize traceback and matrix */

    matrix[0][0] = 0.0;
    traceback[0][0] = 1;            

    for(i=1;i<=lA;i++)
    {
        matrix[0][i] = 0.0; 
        traceback[0][i] = 2;
    }
    for(i=1;i<=lB;i++)
    {
        matrix[i][0] =  0.0;
        traceback[i][0] = 3;
    }

    /* iterate over the sequences */
    
    for(i=1;i<=lB;i++)
    {
        for(j=1;j<=lA;j++)
        {
            
            /* determine the minimum sequence */
            if(i<j)
            {
                minimum = i;
            }
            else
            {
                minimum = j;
            }

            old_score = 0.0;
            old_length = 1;

            /* determine the possible diagonals */
            for(k=0;k<minimum;k++)
            {
                new_score = matrix[i-k-1][j-k-1];

                for(l=k;l>=0;l--)
                {
                        new_score = new_score + scorer[i-l][j-l];
                }
                new_length = k+1;

                if(new_score > old_score)
                {
                    old_score = new_score;
                    old_length = new_length;
                }
            }

            /* determine the previous scores for broken diagonals */
            
            if(arrayB[i-1] < 0 && arrayA[j-1] > 0 && j != lA)
            {
                scoreA = matrix[i-1][j] - 1000000;
            }
            else
            {
                scoreA = matrix[i-1][j];
            }
            if(arrayA[j-1] < 0 && arrayB[i-1] > 0 && i != lB)
            {
                scoreB = matrix[i][j-1] - 1000000;
            }
            else
            {
                scoreB = matrix[i][j-1];
            }

            /* determine the maximum score */
            if(scoreA >= old_score && scoreA > scoreB)
            {
                max_score = scoreA;
                traceback[i][j] = 3;
            }
            else if(old_score > scoreB)
            {
                max_score = old_score;
                for(k=0;k<old_length;k++)
                {
                    traceback[i-k][j-k] = 1;
                }
            }
            else
            {
                max_score = scoreB;
                traceback[i][j] = 2;
            }

            matrix[i][j] = max_score;
        }
    }

    /* get the score for the alignment */
    
    sim = matrix[lB][lA];

    /* carry out the traceback */
    
    i--;
    j--;
    
    while(i > 0 || j > 0)
    {
        if(traceback[i][j] == 3)
        {
            almA[j] = almA[j] + 1;
            i--;
        }
        else if(traceback[i][j] == 1)
        {
            i--;
            j--;
        }
        else
        {
            almB[i] = almB[i] + 1;
            j--;
        }
    }

    return sim;
}

tuple align_pairwise(
        list seqA,
        list seqB,
        list wghA,
        list wghB,
        list resA,
        list resB,
        std::string prsA,
        std::string prsB,
        dict score_dict,
        float scale,
        float sonority_factor,
        std::string mode
        )
{
    const int lA = len(seqA);
    const int lB = len(seqB);

    list outA = extract<list>(seqA.slice(0,lA));
    list outB = extract<list>(seqB.slice(0,lB));

    //std::string vow="Vv<>";
    //bool is_similar;

    int k,l;
    float score,sim;
    
    // define the scoring function
    float (*align)(int *,int,int *,int,float **,float,int *,int *);

    if(mode == "global")
    {
        align = *_global;
    }
    else if(mode == "local")
    {
        align = *_local;
    }
    else if(mode == "repeats")
    {
        align = *_repeats;
    }
    else if(mode == "overlap")
    {
        align = *_overlap;
    }
    else if(mode == "dialign")
    {
        align = *_dialign;
    }

    // create the scorer
    float ** scorer = new float * [lB+1];
    for(k=0;k<lB+1;k++)
    {
        scorer[k] = new float [lA+1];
    }
    scorer[0][0] = 0.0;

    // fill in the scorer with weights for the gaps
    for(k=1;k<lB+1;k++)
    {
        scorer[k][0] = extract<float>(wghB[k-1]);
    }
    for(k=1;k<lA+1;k++)
    {
        scorer[0][k] = extract<float>(wghA[k-1]);
    }

    // fill in the scorer with scores for matchings
    for(k=1;k<lB+1;k++)
    {
        for(l=1;l<lA+1;l++)
        {
           score = extract<float>(score_dict[make_tuple(seqA[l-1],seqB[k-1])]);

	   //if(vow.find(prsA[l-1]) != std::string::npos && vow.find(prsB[k-1]) != std::string::npos)
	   //{
	   //        is_similar = true;
	   //}
	   //else if(vow.find(prsA[l-1]) == std::string::npos && vow.find(prsB[k-1]) == std::string::npos)
	   //{
	   //        is_similar = true;
	   //}
	   //else
	   //{
	   //        is_similar = false;
	   //}
	
           if(prsA[l-1] == prsB[k-1])
           {
                   score = score * (1.0 + sonority_factor);
           }
           scorer[k][l] = score;
        }
    }
    
    // create the arrays passed to _align
    int * arrayA = new int [lA];
    int * arrayB = new int [lB];

    for(k=0;k<lA;k++)
    {
        arrayA[k] = extract<int>(resA[k]);
    }
    for(k=0;k<lB;k++)
    {
        arrayB[k] = extract<int>(resB[k]);
    }
    
    // create almA and almB
    int * almA = new int [lA+1];
    int * almB = new int [lB+1];

    for(k=0;k<lA+1;k++)
    {
        almA[k] = 0;
    }
    for(k=0;k<lB+1;k++)
    {
        almB[k] = 0;
    }
    
    // carry out the alignment
    sim = align(
            arrayA,
            lA,
            arrayB,
            lB,
            scorer,
            scale,
            almA,
            almB
            );
    // retrieve the alignment

    for(k=lA;k>=0;k--)
    {
        if(almA[k] > 0)
        {
            for(l=0;l<almA[k];l++)
            {
                outA.insert(k,"-");
            }
        }
        else if(almA[k] < 0)
        {
            outA[k] = "*";
        }
    }

    for(k=lB;k>=0;k--)
    {
        if(almB[k] > 0)
        {
            for(l=0;l<almB[k];l++)
            {
                outB.insert(k,"-");
            }
        }
        else if(almB[k] < 0)
        {
            outB[k] = '*';
        }
    }

    tuple t = make_tuple(outA,outB,sim);
    
    return t;
}

list align_sequences_pairwise(
        list seqs,
        list weights,
        list restrictions,
        list prosodics,
        dict score_dict,
        float scale,
        float sonority_factor,
        std::string mode
        )
{
    // determine the length of the sequence list
    const int lS = len(seqs);
    
    int i,j,k,l;
    float score;
    
    // this is needed for the clarification of similarity
    //std::string vow="Vv<>";
    //bool is_similar;

    float sim;

    // create the output list
    list alignments = list();

    // define the scoring function
    float (*align)(int *,int,int *,int,float **,float,int *,int *);

    if(mode == "global")
    {
        align = *_global;
    }
    else if(mode == "local")
    {
        align = *_local;
    }
    else if(mode == "repeats")
    {
        align = *_repeats;
    }
    else if(mode == "overlap")
    {
        align = *_overlap;
    }
    else if(mode == "dialign")
    {
        align = *_dialign;
    }
    
    // start the loop over the lists
    for(i=0;i<lS;i++)
    {
        for(j=0;j<lS;j++)
        {
            if(i<j)
            {
                list seqA = extract<list>(seqs[i]);
                list seqB = extract<list>(seqs[j]);
                list wghA = extract<list>(weights[i]);
                list wghB = extract<list>(weights[j]);
                list resA = extract<list>(restrictions[i]);
                list resB = extract<list>(restrictions[j]);
                std::string prsA = extract<std::string>(prosodics[i]);
                std::string prsB = extract<std::string>(prosodics[j]);

                int lA = len(seqA);
                int lB = len(seqB);

                list outA = extract<list>(seqA.slice(0,lA));
                list outB = extract<list>(seqB.slice(0,lB));
                
                // create the scorer
                float ** scorer = new float * [lB+1];
                for(k=0;k<lB+1;k++)
                {
                    scorer[k] = new float [lA+1];
                }
                scorer[0][0] = 0.0;

                // fill in the scorer with weights for the gaps
                for(k=1;k<lB+1;k++)
                {
                    scorer[k][0] = extract<float>(wghB[k-1]);
                }
                for(k=1;k<lA+1;k++)
                {
                    scorer[0][k] = extract<float>(wghA[k-1]);
                }

                // fill in the scorer with scores for matchings
                for(k=1;k<lB+1;k++)
                {
                    for(l=1;l<lA+1;l++)
                    {
                        score = extract<float>(score_dict[make_tuple(seqB[k-1],seqA[l-1])]);
	   		//if(vow.find(prsA[l-1]) != std::string::npos && vow.find(prsB[k-1]) != std::string::npos)
	   		//{
	   		//        is_similar = true;
	   		//}
	   		//else if(vow.find(prsA[l-1]) == std::string::npos && vow.find(prsB[k-1]) == std::string::npos)
	   		//{
	   		//        is_similar = true;
	   		//}
	   		//else
	   		//{
	   		//        is_similar = false;
	   		//}
	
           		//if(prsA[l-1] != prsB[k-1] && is_similar == true)
           		//{
           		//        score = score * (1.0 - sonority_factor);
           		//}
	
            		//if(prsA[l-1] != prsB[k-1] && is_similar == true)
	    		//    	    
	    		//    	    //std::tolower(prsA[l-1]) == std::tolower(prsB[k-1]))
            		//{
            		//        score = score * (1.0 - sonority_factor);
            		//}

           		if(prsA[l-1] == prsB[k-1])
           		{
           		        score = score * (1.0 + sonority_factor);
           		}

                        scorer[k][l] = score;
                    }
                }
                
                // create the arrays passed to _align
                int * arrayA = new int [lA];
                int * arrayB = new int [lB];

                for(k=0;k<lA;k++)
                {
                    arrayA[k] = extract<int>(resA[k]);
                }
                for(k=0;k<lB;k++)
                {
                    arrayB[k] = extract<int>(resB[k]);
                }
                
                // create almA and almB
                int * almA = new int [lA+1];
                int * almB = new int [lB+1];

                for(k=0;k<lA+1;k++)
                {
                    almA[k] = 0;
                }
                for(k=0;k<lB+1;k++)
                {
                    almB[k] = 0;
                }
                
                // carry out the alignment
                sim = align(
                        arrayA,
                        lA,
                        arrayB,
                        lB,
                        scorer,
                        scale,
                        almA,
                        almB
                        );
                // retrieve the alignment

                for(k=lA;k>=0;k--)
                {
                    if(almA[k] > 0)
                    {
                        for(l=0;l<almA[k];l++)
                        {
                            outA.insert(k,"-");
                        }
                    }
                    else if(almA[k] < 0)
                    {
                        outA[k] = "*";
                    }
                }

                for(k=lB;k>=0;k--)
                {
                    if(almB[k] > 0)
                    {
                        for(l=0;l<almB[k];l++)
                        {
                            outB.insert(k,"-");
                        }
                    }
                    else if(almB[k] < 0)
                    {
                        outB[k] = '*';
                    }
                }
                
                alignments.append(make_tuple(outA,outB,sim));

                delete[] arrayA;
                delete[] arrayB;
                for(k=0;k<lB+1;k++)
                {
                    delete[] scorer[k];
                }
                delete[] scorer;
                delete[] almA;
                delete[] almB;
            }
        }
    }
    return alignments;
}

list align_sequence_pairs(
        list seqs,
        list weights,
        list restrictions,
        list prosodics,
        dict score_dict,
        float scale,
        float sonority_factor,
        std::string mode
        )
{
    // determine the length of the sequence list
    const int lS = len(seqs);
    
    int i,k,l; //@check j
    float score;

    float sim;
    
    // specific stuff for the assertion of similarity among prosodic strings
    //std::string vow="Vv<>";
    //bool is_similar;
    
    // create the output list
    list alignments = list();

    // define the scoring function
    float (*align)(int *,int,int *,int,float **,float,int *,int *);

    if(mode == "global")
    {
        align = *_global;
    }
    else if(mode == "local")
    {
        align = *_local;
    }
    else if(mode == "repeats")
    {
        align = *_repeats;
    }
    else if(mode == "overlap")
    {
        align = *_overlap;
    }
    else if(mode == "dialign")
    {
        align = *_dialign;
    }
    
    // start the loop over the lists
    for(i=0;i<lS;i++)
    {
        list seqA = extract<list>(seqs[i][0]);
        list seqB = extract<list>(seqs[i][1]);
        list wghA = extract<list>(weights[i][0]);
        list wghB = extract<list>(weights[i][1]);
        list resA = extract<list>(restrictions[i][0]);
        list resB = extract<list>(restrictions[i][1]);
        std::string prsA = extract<std::string>(prosodics[i][0]);
        std::string prsB = extract<std::string>(prosodics[i][1]);
        
        int lA = len(seqA);
        int lB = len(seqB);

        list outA = extract<list>(seqA.slice(0,lA));
        list outB = extract<list>(seqB.slice(0,lB));
        
        // create the scorer
        float ** scorer = new float * [lB+1];
        for(k=0;k<lB+1;k++)
        {
            scorer[k] = new float [lA+1];
        }
        scorer[0][0] = 0.0;
        
        // fill in the scorer with weights for the gaps
        for(k=1;k<lB+1;k++)
        {
            scorer[k][0] = extract<float>(wghB[k-1]);
        }
        for(k=1;k<lA+1;k++)
        {
            scorer[0][k] = extract<float>(wghA[k-1]);
        }
        // fill in the scorer with scores for matchings
        for(k=1;k<lB+1;k++)
        {
            for(l=1;l<lA+1;l++)
            {
                score = extract<float>(score_dict[make_tuple(seqA[l-1],seqB[k-1])]);
	   	//if(vow.find(prsA[l-1]) != std::string::npos && vow.find(prsB[k-1]) != std::string::npos)
	   	//{
	   	//        is_similar = true;
	   	//}
	   	//else if(vow.find(prsA[l-1]) == std::string::npos && vow.find(prsB[k-1]) == std::string::npos)
	   	//{
	   	//        is_similar = true;
	   	//}
	   	//else
	   	//{
	   	//        is_similar = false;
	   	//}
	
           	//if(prsA[l-1] != prsB[k-1] && is_similar == true)
           	//{
           	//        score = score * (1.0 - sonority_factor);
           	//}
	
            	//if(prsA[l-1] != prsB[k-1] && is_similar == true)
	    	//    	    
	    	//    	    //std::tolower(prsA[l-1]) == std::tolower(prsB[k-1]))
            	//{
            	//        score = score * (1.0 - sonority_factor);
            	//}
           	if(prsA[l-1] == prsB[k-1])
           	{
           	        score = score * (1.0 + sonority_factor);
           	}
                scorer[k][l] = score;
            }
        }
        
        // create the arrays passed to _align
        int * arrayA = new int [lA];
        int * arrayB = new int [lB];

        for(k=0;k<lA;k++)
        {
            arrayA[k] = extract<int>(resA[k]);
        }
        for(k=0;k<lB;k++)
        {
            arrayB[k] = extract<int>(resB[k]);
        }
        
        // create almA and almB
        int * almA = new int [lA+1];
        int * almB = new int [lB+1];

        for(k=0;k<lA+1;k++)
        {
            almA[k] = 0;
        }
        for(k=0;k<lB+1;k++)
        {
            almB[k] = 0;
        }
        
        // carry out the alignment
        sim = align(
                arrayA,
                lA,
                arrayB,
                lB,
                scorer,
                scale,
                almA,
                almB
                );
        // retrieve the alignment

        for(k=lA;k>=0;k--)
        {
            if(almA[k] > 0)
            {
                for(l=0;l<almA[k];l++)
                {
                    outA.insert(k,"-");
                }
            }
            else if(almA[k] < 0)
            {
                outA[k] = "*";
            }
        }

        for(k=lB;k>=0;k--)
        {
            if(almB[k] > 0)
            {
                for(l=0;l<almB[k];l++)
                {
                    outB.insert(k,"-");
                }
            }
            else if(almB[k] < 0)
            {
                outB[k] = '*';
            }
        }
        
        alignments.append(make_tuple(outA,outB,sim));

        delete[] arrayA;
        delete[] arrayB;
        for(k=0;k<lB+1;k++)
        {
            delete[] scorer[k];
        }
        delete[] scorer;
        delete[] almA;
        delete[] almB;
        
    }
    return alignments;
}

dict random_align_sequence_pairs(
        list seqs,
        list weights,
        list restrictions,
        list prosodics,
        dict score_dict,
        float scale,
        float sonority_factor,
        std::string mode,
        int runs
        )
{
    // determine the length of the sequence list
    const int lS = len(seqs);
    
    int i,k,l,m,n,lA,lB; //@check j
    float score;

    float sim;

    list outA,outB,seqA,seqB,wghA,wghB,resA,resB;
    std::string prsA,prsB;

    // create the output dictionary
    dict corrs = dict();
    tuple tmp_pair;
    
    // specific stuff for the assertion of similarity among prosodic strings
    //std::string vow="Vv<>";
    //bool is_similar;

    // define the scoring function
    float (*align)(int *,int,int *,int,float **,float,int *,int *);

    if(mode == "global")
    {
        align = *_global;
    }
    else if(mode == "local")
    {
        align = *_local;
    }
    else if(mode == "repeats")
    {
        align = *_repeats;
    }
    else if(mode == "overlap")
    {
        align = *_overlap;
    }
    else if(mode == "dialign")
    {
        align = *_dialign;
    }
    
    // create the vector of numbers by which the sequences are shuffled
    std::vector<int> v;

    for (int i=0;i<lS;i++)
        {
            v.push_back(i);
        }
    
    // create a dictionary which stores all aligned pairs in order to avoid
    // that alignments, carried out once, are repeated later on
    dict alm_pairs = dict();

    // start a loop according to the number of runs
    for(m=0;m<runs;m++)
    {
        std::random_shuffle(v.begin(),v.end());
    
        // start the loop over the lists
        for(i=0;i<lS;i++)
        {
            // check whether the alignment has been carried out before
            if(alm_pairs.has_key(make_tuple(v[i],i)) != true)
            {
                // carry out the regular alignment 
                seqA = extract<list>(seqs[v[i]][0]);
                seqB = extract<list>(seqs[i][1]);
                wghA = extract<list>(weights[v[i]][0]);
                wghB = extract<list>(weights[i][1]);
                resA = extract<list>(restrictions[v[i]][0]);
                resB = extract<list>(restrictions[i][1]);
                prsA = extract<std::string>(prosodics[v[i]][0]);
                prsB = extract<std::string>(prosodics[i][1]);
                
                lA = len(seqA);
                lB = len(seqB);

                outA = extract<list>(seqA.slice(0,lA));
                outB = extract<list>(seqB.slice(0,lB));
                

                // create the scorer
                float ** scorer = new float * [lB+1];
                for(k=0;k<lB+1;k++)
                {
                    scorer[k] = new float [lA+1];
                }
                scorer[0][0] = 0.0;
                
                // fill in the scorer with weights for the gaps
                for(k=1;k<lB+1;k++)
                {
                    scorer[k][0] = extract<float>(wghB[k-1]);
                }
                for(k=1;k<lA+1;k++)
                {
                    scorer[0][k] = extract<float>(wghA[k-1]);
                }
                // fill in the scorer with scores for matchings
                for(k=1;k<lB+1;k++)
                {
                    for(l=1;l<lA+1;l++)
                    {
                        score = extract<float>(score_dict[make_tuple(seqA[l-1],seqB[k-1])]);
	   		//if(vow.find(prsA[l-1]) != std::string::npos && vow.find(prsB[k-1]) != std::string::npos)
	   		//{
	   		//        is_similar = true;
	   		//}
	   		//else if(vow.find(prsA[l-1]) == std::string::npos && vow.find(prsB[k-1]) == std::string::npos)
	   		//{
	   		//        is_similar = true;
	   		//}
	   		//else
	   		//{
	   		//        is_similar = false;
	   		//}
	
           		//if(prsA[l-1] != prsB[k-1] && is_similar == true)
           		//{
           		//        score = score * (1.0 - sonority_factor);
           		//}
	
           		if(prsA[l-1] == prsB[k-1])
           		{
           		        score = score * (1.0 + sonority_factor);
           		}
                        scorer[k][l] = score;
                    }
                }
                
                // create the arrays passed to _align
                int * arrayA = new int [lA];
                int * arrayB = new int [lB];

                for(k=0;k<lA;k++)
                {
                    arrayA[k] = extract<int>(resA[k]);
                }
                for(k=0;k<lB;k++)
                {
                    arrayB[k] = extract<int>(resB[k]);
                }
                
                // create almA and almB
                int * almA = new int [lA+1];
                int * almB = new int [lB+1];

                for(k=0;k<lA+1;k++)
                {
                    almA[k] = 0;
                }
                for(k=0;k<lB+1;k++)
                {
                    almB[k] = 0;
                }
                
                // carry out the alignment
                sim = align(
                        arrayA,
                        lA,
                        arrayB,
                        lB,
                        scorer,
                        scale,
                        almA,
                        almB
                        );
                // retrieve the alignment

                for(k=lA;k>=0;k--)
                {
                    if(almA[k] > 0)
                    {
                        for(l=0;l<almA[k];l++)
                        {
                            outA.insert(k,"-");
                        }
                    }
                    else if(almA[k] < 0)
                    {
                        outA[k] = "*";
                    }
                }

                for(k=lB;k>=0;k--)
                {
                    if(almB[k] > 0)
                    {
                        for(l=0;l<almB[k];l++)
                        {
                            outB.insert(k,"-");
                        }
                    }
                    else if(almB[k] < 0)
                    {
                        outB[k] = '*';
                    }
                }
                // delete the arrays of the alignment analysis
                delete[] arrayA;
                delete[] arrayB;
                for(k=0;k<lB+1;k++)
                {
                    delete[] scorer[k];
                }
                delete[] scorer;
                delete[] almA;
                delete[] almB;

                // introduce the alignment into the dictionary which stores
                // everything
                alm_pairs[make_tuple(v[i],i)] = make_tuple(outA,outB);
            }
            else
            {
                // get the respective alignment
                outA = extract<list>(alm_pairs[make_tuple(v[i],i)][0]);
                outB = extract<list>(alm_pairs[make_tuple(v[i],i)][1]);
            }
            // fill in the values for the correspondence dictionary, mind to
            // make a specific procedure when filling in local alignments!
            if(mode != "local")
            {

                for(n=0;n<len(outA);n++)
                {
                    tmp_pair = make_tuple(outA[n],outB[n]);

                    if(corrs.has_key(tmp_pair) == false)
                    {
                        corrs[tmp_pair] = 1.0 / runs;
                    }
                    else
                    {
                        corrs[tmp_pair] += 1.0 / runs;
                    }

                }
            }
            else
            {
                list localA = list();
                list localB = list();

                for(n=0;n<len(outA);n++)
                {
                    if(outA[n] != "*")
                    {
                        localA.append(outA[n]);
                    }
                }
                for(n=0;n<len(outB);n++)
                {
                    if(outB[n] != "*")
                    {
                        localB.append(outB[n]);
                    }
                }
                for(n=0;n<len(localA);n++)
                {
                    tmp_pair = make_tuple(localA[n],localB[n]);

                    if(corrs.has_key(tmp_pair) == false)
                    {
                        corrs[tmp_pair] = 1.0 / runs;
                    }
                    else
                    {
                        corrs[tmp_pair] += 1.0 / runs;
                    }
                }
            }
            
        }
    }
    return corrs;
}

/* Function returns the normalized Levenshtein-distance between two lists. 
 */
float edit_dist(
        list seqA,
        list seqB
        )
{
    
    // set initial values
    const int m = len(seqA);
    const int n = len(seqB);

    // create the matrix
    int ** matrix = new int *[n+1];
    for(int i=0;i<=n;i++)
    {
        matrix[i] = new int [m+1];
    }
    
    // initialize the matrix
    matrix[0][0] = 0;
    int i=0, j=0;
    for(j=1;j<=m;j++)
    {
        matrix[0][j] = j;
    }
    for(i=1;i<=n;i++)
    {
        matrix[i][0] = i;
    }

    // set the parameters
    int penalty;
    int gapA, gapB, match;
    int min;
    float dist;
    
    // fill in the matrix
    for(i=1;i<=n;i++)
    {
        for(j=1;j<=m;j++)
        {
            if(seqA[j-1] == seqB[i-1])
            {
                penalty = 0;
            }
            else
            {
                penalty = 1;
            }
            gapA = matrix[i-1][j] + 1;
            match = matrix[i-1][j-1] + penalty;
            gapB = matrix[i][j-1] + 1;

            // get the minimum value
            if(gapA <= match && gapA <= gapB)
            {
                min = gapA;
            }
            else if(match < gapB)
            {
                min = match;
            }
            else
            {
                min = gapB;
            }

            matrix[i][j] = min;
        }
    }    
    dist = (matrix[n][m] * 1.0) / (std::max(n,m) * 1.0);
    
    return dist;
}


BOOST_PYTHON_MODULE(align)
{
    def(
            "align_sequences_pairwise",
            &align_sequences_pairwise
       );
    def(
            "align_pairwise",
            &align_pairwise
       );
    def(
            "align_sequence_pairs",
            &align_sequence_pairs
       );
    def(
            "random_align_sequence_pairs",
            &random_align_sequence_pairs
       );
    def(
            "edit_dist",
            &edit_dist
       );
}
