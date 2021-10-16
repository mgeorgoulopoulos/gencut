#ifndef __GENE_COLLECTIONS_H__
#define __GENE_COLLECTIONS_H__

#include <utils/AtomBox.h>

#include <set>
#include <vector>

using GeneId = Atom;

class GeneSet : public std::set<GeneId> {};

class GeneList : public std::vector<GeneId> {};

#endif // __GENE_COLLECTIONS_H__