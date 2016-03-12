from __future__ import division, unicode_literals, print_function
import re

from lingpy.thirdparty import LoadTree as Tree


class _TreeDist(object):
    """
    Private class with static methods for the computation of RF-distance and
    Normalized RF distance.
    """

    # basic constants
    LBRA = re.compile("\(")
    RBRA = re.compile("\)")

    @staticmethod
    def grf(treeA, treeB, distance='grf'):
        """
        Computes the generalized Robinson fould distance between two trees.
        """

        # prepare the trees [probably not necessary @lingulist]
        for tree in [treeA, treeB]:
            for old, new in [(";", ""), ("/", "-"), ("\n", "")]:
                tree = tree.replace(old, new)
        treeA = treeA.replace(" ", "_")

        # get lingpy-trees from treeA and treeB
        ntreeAtree = Tree(treeA + ';')

        lang_settreeA = set(ntreeAtree.taxa)

        ntreeBtree = Tree(treeB + ';')
        lang_settreeB = set(ntreeBtree.taxa)

        # check for identical number of taxa
        if len(lang_settreeA) != len(lang_settreeB):
            raise ValueError(
                "The number of taxonomic units should be identical in both trees!")

        treeA_parts, l_treeA = _TreeDist.get_bipartition(treeA)
        treeB_parts, l_etn = _TreeDist.get_bipartition(treeB)

        # calculate partition distance (= symmetric difference) using
        # lingpy tree class)
        e = 0.0
        i_treeA = len(treeA_parts)
        i_treeB = len(treeB_parts)

        for upart in treeA_parts.keys():
            upart1 = l_treeA - upart
            if upart in treeB_parts or upart1 in treeB_parts:
                e += 1.0

        e_mod = 0.0
        f = 0
        for upart in treeA_parts.keys():
            upart1 = l_treeA - upart
            emod = None
            f += 1
            for epart in treeB_parts.keys():
                if upart <= epart or upart1 <= epart:
                    emod = True
                else:
                    epart1 = l_etn - epart
                    if upart <= epart1 or upart1 <= epart1:
                            emod = True
                    else:
                            emod = False
                            break
                if emod == False:
                    break
            if emod:
                e_mod += 1.0

        grf = (i_treeA - e_mod) / i_treeA
        rf = (i_treeA + i_treeB - 2 * e) / (i_treeA + i_treeB)
        return grf if distance == 'grf' else rf

    @staticmethod
    def get_bipartition(tree):
        partition_list = []
        temp_stack = []
        ind_list = []
        hash_lang = {}  # nltk.defaultdict(int)
        lang_cnt = 0
        tree_list = tree.split(",")
        for i, elem in enumerate(tree_list):
            if elem.find("(") > -1:
                for k in [m.start() for m in re.finditer(_TreeDist.LBRA, elem)]:
                    ind_list.append(i)
                lang = elem.strip("(")
                lang = lang.split(":")[0]
                lang = lang.replace("(", "-")
                lang = lang.replace(")", "-")
                lang = lang.replace(" ", "_").replace("'", "")
                if lang in hash_lang:
                    temp_stack.append(hash_lang[lang])
                else:
                    lang_cnt += 1
                    hash_lang[lang] = str(lang_cnt)
                    hash_lang[lang] = lang
                    temp_stack.append(hash_lang[lang])
            elif elem.find(")") > -1:
                lang = elem.replace(")", "")
                lang = lang.split(":")[0].strip()
                if lang in hash_lang:
                    temp_stack.append(hash_lang[lang])
                else:
                    lang_cnt += 1
                    hash_lang[lang] = str(lang_cnt)
                    hash_lang[lang] = lang
                    temp_stack.append(hash_lang[lang])
                for k in [m.start() for m in re.finditer(_TreeDist.RBRA, elem)]:
                    partition = temp_stack[ind_list.pop():]
                    p1 = partition
                    partition_list.append(p1)
            else:
                lang = elem.split(":")[0]
                if lang in hash_lang:
                    temp_stack.append(hash_lang[lang])
                else:
                    lang_cnt += 1
                    hash_lang[lang] = str(lang_cnt)
                    hash_lang[lang] = lang
                    temp_stack.append(hash_lang[lang])

        if len(ind_list) > 0:
            raise ValueError("Cannot compute the bipartition!")

        lang_set = frozenset(partition_list[-1])
        final_parts = {}
        for x in partition_list:
            set_x = frozenset(x)
            set_x1 = lang_set - set_x
            if len(set_x1) == 1 or len(set_x) == 1:
                continue

            if len(set_x) > 0 and len(set_x1) > 0:
                if set_x not in final_parts and set_x1 not in final_parts:
                    final_parts[set_x] = True

        return final_parts, lang_set
