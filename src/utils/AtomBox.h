#ifndef __ATOM_BOX_H__
#define __ATOM_BOX_H__

#include <map>
#include <string>

using Atom = int;

class AtomBox : public std::map<std::string, Atom> {
  public:
	bool contains(const std::string &geneName) const {
		return find(geneName) != end();
	}

	Atom atom(const std::string &geneName) {
		auto it = find(geneName);
		if (it != end())
			return it->second;

		// Not found - insert
		const Atom atom = (int)size() + 1;
		(*this)[geneName] = atom;

		return atom;
	}

	// Too bored to do this properly now
	std::string name(Atom atom) {
		for (auto it = begin(); it != end(); it++) {
			if (it->second == atom)
				return it->first;
		}
		return "";
	}
};

#endif // __ATOM_BOX_H__