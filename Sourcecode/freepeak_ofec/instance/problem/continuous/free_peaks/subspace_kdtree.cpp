#include "subspace_kdtree.h"
#include "param_creator.h"
#include <queue>

namespace ofec::free_peaks {
	void SubspaceKDTree::createTree(const std::vector<std::pair<double, double>>& ranges,
		const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>& tree_name)
	{
		root_box.reset(new Box(tree_name.front().first));
		name_box_subproblem[root_box->name] = make_pair(root_box, nullptr);
		std::shared_ptr<Box> tmp_box;
		for (auto& it : tree_name) {
			for (auto& it2 : it.second) {
				tmp_box.reset(new Box(it2.first, it2.second));
				name_box_subproblem[it2.first] = make_pair(tmp_box, nullptr);
				name_box_subproblem[it.first].first->children.push_back(tmp_box);
			}
		}
		tree.reset(new nanoflann::KDTreeSpace<double>(tree_name.front().second, ranges));
		tree->buildIndex();
		if (tree_name.size() > 1) {
			for (unsigned i(1); i < tree_name.size(); ++i) {
				tree->addSubtree(tree_name[i].first, tree_name[i].second);
			}
		}
	}

	std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> SubspaceKDTree::treeData() const {
		std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>> tmp_tree_data;
		std::queue<std::shared_ptr<Box>> qb;
		qb.push(root_box);
		std::shared_ptr<Box> tmp;
		std::pair<std::string, std::vector<std::pair<std::string, double>>> tmp_pair;
		while (!qb.empty()) {
			tmp = qb.front();
			qb.pop();
			tmp_pair.first = tmp->name;
			tmp_pair.second.clear();
			for (auto& it : tmp->children) {
				tmp_pair.second.emplace_back(make_pair(it->name, it->ratio));
				if (!it->children.empty()) {
					qb.push(it);
				}
			}
			tmp_tree_data.emplace_back(tmp_pair);
		}
		return tmp_tree_data;
	}

	void SubspaceKDTree::clear() {
		root_box.reset();
		name_box_subproblem.clear();
		tree.reset();
	}
}