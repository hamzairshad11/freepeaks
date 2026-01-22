//*****************************************************************************************************************************************

//part of the source code is from nanoflann lib at https://github.com/jlblancoc/nanoflann

#ifndef  NANOFLANN_HPP_
#define  NANOFLANN_HPP_

#include <vector>
#include <list>
#include <set>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <utility>
#include <iostream>
#include "../../core/definition.h"
#include <map>
#include <memory>
#include "pool_allocator.h"
#include <string>

namespace nanoflann
{
	template <typename ElementType, typename IndexType = size_t>
	class KDTreeSpace
	{
	private:
		/** Hidden copy constructor, to disallow copying indices (Not implemented) */
		KDTreeSpace(const KDTreeSpace<ElementType, IndexType> &) = delete;
	protected:

		/**
		*  Array of indices to vectors in the dataset.
		*/
		std::vector<IndexType> m_vind;

		/**
		* The dataset used by this index
		*/
		std::vector<std::vector<ElementType>> m_point_data; //!< The source of our data
		std::vector<ElementType> m_ratio_data;
		std::vector<std::pair<std::string, ElementType>> m_name_ratio_data;

		size_t m_size; //!< Number of current points in the dataset
		size_t m_dim;  //!< Dimensionality of each data point
		size_t m_relative_depth = 0;

		/*--------------------- Internal Data Structures --------------------------*/
		struct BoundingBox {
			std::vector<std::pair<ElementType, ElementType>> box;
			ofec::Real rat = 1.0, volume;
			size_t depth = 1;
			std::string name;
		};

		struct Node
		{
			/** Union used because a node can be either a LEAF node or a non-leaf node, so both data fields are never used simultaneously */
			union {
				struct {
					IndexType    idx_region;
					BoundingBox *region;
				}lr;
				struct {
					int          divfeat; //!< Dimension used for subdivision.
					ElementType pivot; // pivot value for division
					IndexType    idx_sample;		// index
					ElementType low, high;	//boundary for the box to be cutted.
				}sub;
			};
			size_t depth;
			Node *child1, *child2, *parent;  //!< Child nodes (both=nullptr mean its a leaf node)
			std::string name;
		};
		typedef Node* NodePtr;

		/** The KD-tree used to find regions */
		NodePtr m_root;
		BoundingBox m_init_region;
		std::vector<std::unique_ptr<BoundingBox>> m_subregions;
		std::set<IndexType> m_shelved_indices;
		std::list<NodePtr> m_shelved_nodes;
		std::map<std::string, NodePtr> m_name_to_leaf_node;

		/**
		* Pooled memory allocator.
		*
		* Using a pooled memory allocator is more efficient
		* than allocating memory directly when there is a large
		* number small of memory allocations.
		*/
		PooledAllocator<Node> m_pool;
		ofec::Real m_lrat = 0, m_srat = 1;
		int m_lbox = 0, m_sbox = 0;
		int m_mode;

	public:
		// constructor based on empty data 
		KDTreeSpace() : m_root(nullptr), m_size(0), m_dim(0), m_mode(0) {}

		// constructor based on points
		KDTreeSpace(const std::vector<std::vector<ElementType>> &points,
			const std::vector<std::pair<ElementType, ElementType>> &init_box) :
			m_point_data(points), m_root(nullptr), 
			m_size(points.size()), m_dim(init_box.size()), m_mode(1)
		{
			m_init_region.box = init_box;
			initVindPoint();
		}

		KDTreeSpace(const std::vector<std::vector<ElementType>> &points) :
			m_point_data(points), m_root(nullptr), 
			m_size(points.size()), m_dim(0), m_mode(1)
		{
			initVindPoint();
		}

		// constructor based on ratios
		KDTreeSpace(const std::vector<ElementType> &ratios,
			const std::vector<std::pair<ElementType, ElementType>> &init_box) :
			m_ratio_data(ratios), m_root(nullptr), 
			m_size(ratios.size()), m_dim(init_box.size()), m_mode(2)
		{
			m_init_region.box = init_box;
			initVindRatio();
		}

		KDTreeSpace(const std::vector<ElementType> &ratios) :
			m_ratio_data(ratios), m_root(nullptr), 
			m_size(ratios.size()), m_dim(0), m_mode(2)
		{
			initVindRatio();
		}

		//constructor based on names and ratios
		KDTreeSpace(const std::vector<std::pair<std::string, ElementType>> &names_ratios,
			const std::vector<std::pair<ElementType, ElementType>> &init_box) :
			m_name_ratio_data(names_ratios), m_root(nullptr),
			m_size(names_ratios.size()), m_dim(init_box.size()), m_mode(3)
		{
			m_init_region.box = init_box;
			initVindNameRatio();
		}

		KDTreeSpace(const std::vector<std::pair<std::string, ElementType>> &names_ratios) :
			m_name_ratio_data(names_ratios), m_root(nullptr),
			m_size(names_ratios.size()), m_dim(0), m_mode(3)
		{
			initVindNameRatio();
		}

		// change the initial region or data 
		void setInitBox(const std::vector<std::pair<ElementType, ElementType>> &init_box) {
			m_init_region.box = init_box;
			m_dim = init_box.size();
		}

		const std::vector<std::pair<ElementType, ElementType>>& getInitBox() {
			return m_init_region.box;
		}

		// change the point data 
		void inputPointData(const std::vector<std::vector<ElementType>> &points) {
			m_point_data = points;
			m_size = points.size();
			m_mode = 1;
			initVindPoint();
		}

		// change the ratio data 
		void inputRatioData(const std::vector<ElementType> &ratios) {
			m_ratio_data = ratios;
			m_size = ratios.size();
			m_mode = 2;
			initVindRatio();
		}

		// change the name and ratio data 
		void inputNameRatioData(const std::vector<std::pair<std::string, ElementType>> &names_ratios) {
			m_name_ratio_data = names_ratios;
			m_size = names_ratios.size();
			m_mode = 3;
			initVindRatio();
		}

		std::vector<IndexType> regionsIndices() const {
			std::vector<IndexType> idxs;
			for (auto i = 0; i < m_subregions.size(); ++i)
				if (m_shelved_indices.count(i) == 0)
					idxs.push_back(i);
			return idxs;
		}

		//Printout the idxBoxes
		void regionShow() {
			for (auto i = 0; i < m_subregions.size(); ++i) {
				if (m_shelved_indices.count(i) == 0) {
					std::cout << i << ":   ";
					for (auto j = 0; j < m_dim; ++j) {
						std::cout << "(" << m_subregions[i]->box[j].first << " , " << m_subregions[i]->box[j].second << ") ";
					}
					std::cout << std::endl;
				}
			}
		}

		/** Standard destructor */
		~KDTreeSpace() {}

		/** Frees the previously-built index. Automatically called within buildIndex(). */
		void freeIndex() {
			m_pool.free_all();
			m_root = nullptr;
			m_subregions.clear();
			m_shelved_indices.clear();
			m_shelved_nodes.clear();
			m_name_to_leaf_node.clear();
		}

		/**
		* Builds the index
		*/
		void buildIndex() {
			freeIndex();
			if (m_mode == 1) { // randomly constr
				initVindPoint();
				m_root = divideByPoint(0, m_size, m_init_region);   // construct the tree
			}
			else if (m_mode == 2) {
				initVindRatio();
				m_root = divideByRatio(0, m_size, m_init_region);
			}
			else if (m_mode == 3) {
				initVindNameRatio();
				m_root = divideByNameRatio(0, m_size, m_init_region);
			}
		}

		/** Returns number of leaf nodes  */
		size_t size() const {
			if (m_mode == 1) return m_size + 1;
			else return m_size;
		}

		size_t veclen() const { return m_dim; }
		int smallestBox() const { return m_sbox; }
		int largestBox() const { return m_lbox; }
		void setRelativeDepth(size_t depth) { m_relative_depth = depth; }

		const std::vector<std::pair<ElementType, ElementType>>& getBox(IndexType idx) const { 
			return m_subregions.at(idx)->box; 
		}

		const std::vector<std::pair<ElementType, ElementType>>& getBox(const std::string name) const { 
			return m_name_to_leaf_node.at(name)->lr.region->box;
		}

		ofec::Real getBoxVolume(IndexType idx) const { 
			return m_subregions.at(idx)->volume; 
		}

		size_t getDepth(IndexType idx) const { 
			return m_subregions.at(idx)->depth; 
		}

		const std::string& getName(IndexType idx) const {
			return m_subregions.at(idx)->name;
		}

		int getRegionIdx(const std::vector<ElementType> &p) const { 
			return enqury(p, m_root)->lr.idx_region; 
		}

		const std::string& getRegionName(const std::vector<ElementType> &p) const {
			return enqury(p, m_root)->name; 
		}

		int splitRegion(IndexType idx, int *dim = nullptr, ElementType *split_pivot = nullptr) {
			NodePtr node = nullptr;
			leafNode(idx, m_root, node);
			if (node == nullptr) return -1;
			NodePtr node1, node2;
			if (!m_shelved_nodes.empty()) {
				node1 = m_shelved_nodes.front();
				m_shelved_nodes.pop_front();
			}
			else
				node1 = m_pool.allocate();
			if (!m_shelved_nodes.empty()) {
				node2 = m_shelved_nodes.front();
				m_shelved_nodes.pop_front();
			}
			else
				node2 = m_pool.allocate();
			node1->depth = node->depth + 1;
			node2->depth = node->depth + 1;
			node1->parent = node;
			node2->parent = node;
			node1->child1 = node1->child2 = nullptr;
			node2->child1 = node2->child2 = nullptr;
			node->child1 = node1;
			node->child2 = node2;
			if (dim != nullptr)
				node->sub.divfeat = *dim;
			else
				node->sub.divfeat = node->depth % m_dim;
			node->sub.low = m_subregions[idx]->box[node->sub.divfeat].first;
			node->sub.high = m_subregions[idx]->box[node->sub.divfeat].second;
			if (!split_pivot || *split_pivot < node->sub.low || *split_pivot > node->sub.high)
				node->sub.pivot = (node->sub.low + node->sub.high) / 2;
			else
				node->sub.pivot = *split_pivot;
			node1->lr.idx_region = idx;
			m_subregions[node1->lr.idx_region]->depth++;
			if (m_shelved_indices.empty()) {
				node2->lr.idx_region = m_subregions.size();
				m_subregions.emplace_back(new BoundingBox(*m_subregions[node1->lr.idx_region]));
			}
			else {
				node2->lr.idx_region = *m_shelved_indices.begin();
				m_shelved_indices.erase(m_shelved_indices.begin());
				m_subregions[node2->lr.idx_region]->depth = m_subregions[node1->lr.idx_region]->depth;
			}
			m_subregions[node1->lr.idx_region]->box[node->sub.divfeat].second = node->sub.pivot;
			m_subregions[node2->lr.idx_region]->box[node->sub.divfeat].first = node->sub.pivot;
			node1->lr.region = m_subregions[node1->lr.idx_region].get();
			node2->lr.region = m_subregions[node2->lr.idx_region].get();
			m_subregions[node1->lr.idx_region]->name = node1->name;
			m_subregions[node2->lr.idx_region]->name = node2->name;
			auto pre_vol = m_subregions[node1->lr.idx_region]->volume;
			m_subregions[node1->lr.idx_region]->volume = pre_vol * ((node->sub.pivot - node->sub.low) / (node->sub.high - node->sub.low));
			m_subregions[node2->lr.idx_region]->volume = pre_vol - m_subregions[node1->lr.idx_region]->volume;
			m_size++;
			return node2->lr.idx_region;
		}

		void mergeRegion(IndexType idx, std::list<IndexType> &merged_regs) {      // Merge the region and its sibling region
			NodePtr node = nullptr;
			leafNode(idx, m_root, node);
			if (node == nullptr || node == m_root) return;
			mergeNode(node, merged_regs);
		}

		void unionsAtDepth(size_t depth, std::vector<std::vector<IndexType>> &region_unions) const {
			if (!region_unions.empty()) region_unions.clear();
			findUnionsAtDepth(depth, m_root, region_unions);
		}

		void findNeighbor(IndexType idx, std::list<IndexType> &neighbors) const {
			if (!neighbors.empty()) neighbors.clear();
			neighborCheck(idx, m_root, neighbors);
		}

		bool checkAdjacency(int idx1, int idx2) const {
			bool result = true;
			const auto &box1 = m_subregions[idx1]->box;
			const auto &box2 = m_subregions[idx2]->box;
			for (size_t j = 0; j < m_dim; ++j) {
				if (box1[j].second < box2[j].first || box2[j].second < box1[j].first) {
					result = false;
					break;
				}
			}
			//if (result) {
			//	bool all_adj = true;
			//	for (size_t j = 0; j < m_dim; ++j) {
			//		if (box1[j].second != box2[j].first && box2[j].second != box1[j].first) {
			//			all_adj = false;
			//			break;
			//		}
			//	}
			//	result = !all_adj;
			//}
			return result;
		}

		void addSubtree(IndexType idx, KDTreeSpace &subtree, std::map<IndexType, IndexType> &add_id) {
			NodePtr node = nullptr;
			leafNode(idx, m_root, node);
			if (node == nullptr)	return;
			if (node->lr.region->box != subtree.m_init_region.box)	return;
			auto root = subtree.m_root;
			if (root->child1 == nullptr && root->child2 == nullptr) return;
			m_shelved_indices.insert(idx);
			m_size--;
			copyLeafNodes(node, root, add_id);
		}

		bool addSubtree(IndexType idx, const std::vector<ElementType> &ratios) {
			NodePtr node = nullptr;
			leafNode(idx, m_root, node);
			if (node == nullptr)	return false;
			auto box = node->lr.region->box;
			KDTreeSpace subtree(ratios, box);
			std::map<IndexType, IndexType> add_id;
			m_shelved_indices.insert(node->lr.idx_region);
			m_size--;
			copyLeafNodes(node, subtree.m_root, add_id);
			return true;
		}

		bool addSubtree(const std::string &name, const std::vector<std::pair<std::string, ElementType>> &names_ratios) {
			if (!m_name_to_leaf_node.count(name) || names_ratios.size() < 2)
				return false;
			NodePtr node = m_name_to_leaf_node[name];
			auto box = node->lr.region->box;
			KDTreeSpace subtree(names_ratios, box);
			subtree.setRelativeDepth(node->depth);
			subtree.buildIndex();
			std::map<IndexType, IndexType> add_id;
			m_shelved_indices.insert(node->lr.idx_region);
			m_size--;
			copyLeafNodes(node, subtree.m_root, add_id);
			return true;
		}

	private:
		/** Make sure the auxiliary list \a vind has the same size than the current dataset, and re-generate if size has changed. */
		void initVindPoint() {
			// Create a permutable array of indices to the input vectors.
			m_size = m_point_data.size();
			if (m_vind.size() != m_size) m_vind.resize(m_size);
			size_t k = 0;
			for (auto &i : m_vind) i = k++;
		}

		/** Make sure the auxiliary list \a vind has the same size than the current ratiodata, and re-generate if size has changed. */
		void initVindRatio() {
			// Create a permutable array of indices to the input vectors.
			m_size = m_ratio_data.size();
			if (m_vind.size() != m_size) m_vind.resize(m_size);
			size_t k = 0;
			for (auto &i : m_vind) i = k++;
		}

		/** Make sure the auxiliary list \a vind has the same size than the current ratiodata, and re-generate if size has changed. */
		void initVindNameRatio() {
			// Create a permutable array of indices to the input vectors.
			m_size = m_name_ratio_data.size();
			if (m_vind.size() != m_size) m_vind.resize(m_size);
			size_t k = 0;
			for (auto &i : m_vind) i = k++;
		}

		/// Helper accessor to the dataset points:
		ElementType getDataset(size_t idx, int component) const {
			return m_point_data[idx][component];
		}

		/// Helper accessor to the ratiodata:
		ElementType getRatiodata(size_t idx) const {
			return m_ratio_data[idx];
		}

		/**
		* Create a tree node that subdivides the list of vecs from vind[first]
		* to vind[last].  The routine is called recursively on each sublist.
		*
		* @param left index of the first vector
		* @param right index of the last vector
		*/
		NodePtr divideByPoint(const IndexType left, const IndexType right, BoundingBox &region, size_t depth = 0) {
			NodePtr node = m_pool.allocate(); // allocate memory

			/*a leaf node,create a sub-region. */
			if ((right - left) <= 0) {
				node->child1 = node->child2 = nullptr;    /* Mark as leaf node. */
				node->lr.idx_region = m_subregions.size();
				node->depth = depth;
				m_subregions.emplace_back(new BoundingBox(region));
				boxRatio(*m_subregions.back(), m_subregions.size() - 1);
			}
			else {
				IndexType idx;
				int cutfeat;
				ElementType cutval;
				divisionByPoint(&m_vind[0] + left, right - left, idx, cutfeat, cutval, region, depth);

				node->sub.idx_sample = m_vind[left + idx];
				node->sub.divfeat = cutfeat;
				node->sub.low = region.box[cutfeat].first;
				node->sub.high = region.box[cutfeat].second;
				node->depth = depth;

				BoundingBox left_subregion(region);
				left_subregion.box[cutfeat].second = cutval;
				NodePtr temp = divideByPoint(left, left + idx, left_subregion, depth + 1);
				node->child1 = temp;
				temp->parent = node;

				BoundingBox right_subregion(region);
				right_subregion.box[cutfeat].first = cutval;
				temp = divideByPoint(left + idx + 1, right, right_subregion, depth + 1);
				node->child2 = temp;
				temp->parent = node;

				node->sub.pivot = cutval;
			}
			return node;
		}

		NodePtr divideByRatio(const IndexType left, const IndexType right, BoundingBox &region, size_t depth = 0) {
			NodePtr node = m_pool.allocate(); // allocate memory
			if ((right - left) <= 1) {
				node->child1 = node->child2 = nullptr;    /* Mark as leaf node. */
				node->lr.idx_region = m_subregions.size();
				node->depth = depth;
				m_subregions.emplace_back(new BoundingBox(region));
				boxRatio(*m_subregions.back(), m_subregions.size() - 1);
				node->lr.region = m_subregions.back().get();
			}
			else {
				IndexType idx;
				int cutfeat;
				ElementType cutval;
				divisionByRatio(&m_vind[0] + left, right - left, idx, cutfeat, cutval, region, depth);

				node->sub.idx_sample = m_vind[idx];
				node->sub.divfeat = cutfeat;
				node->sub.low = region.box[cutfeat].first;
				node->sub.high = region.box[cutfeat].second;
				node->depth = depth;

				BoundingBox left_subregion(region);
				left_subregion.box[cutfeat].second = cutval;
				left_subregion.depth = depth + 1;
				NodePtr temp = divideByRatio(left, idx, left_subregion, depth + 1);
				node->child1 = temp;
				temp->parent = node;

				BoundingBox right_subregion(region);
				right_subregion.box[cutfeat].first = cutval;
				right_subregion.depth = depth + 1;
				temp = divideByRatio(idx, right, right_subregion, depth + 1);
				node->child2 = temp;
				temp->parent = node;
				node->sub.pivot = cutval;
			}
			return node;
		}

		NodePtr divideByNameRatio(const IndexType left, const IndexType right, BoundingBox &region, size_t depth = 0) {
			NodePtr node = m_pool.allocate(); // allocate memory
			if ((right - left) <= 1) {
				m_name_to_leaf_node[node->name = m_name_ratio_data[left].first] = node;
				node->child1 = node->child2 = nullptr;    /* Mark as leaf node. */
				node->lr.idx_region = m_subregions.size();
				node->depth = depth;
				m_subregions.emplace_back(new BoundingBox(region));
				boxRatio(*m_subregions.back(), m_subregions.size() - 1);
				node->lr.region = m_subregions.back().get();
				m_subregions.back()->name = node->name;
			}
			else {
				IndexType idx;
				int cutfeat;
				ElementType cutval;
				divisionByNameRatio(&m_vind[0] + left, right - left, idx, cutfeat, cutval, region, depth);

				node->sub.idx_sample = m_vind[idx];
				node->sub.divfeat = cutfeat;
				node->sub.low = region.box[cutfeat].first;
				node->sub.high = region.box[cutfeat].second;
				node->depth = depth;

				BoundingBox left_subregion(region);
				left_subregion.box[cutfeat].second = cutval;
				left_subregion.depth = depth + 1;
				NodePtr temp = divideByNameRatio(left, idx, left_subregion, depth + 1);
				node->child1 = temp;
				temp->parent = node;

				BoundingBox right_subregion(region);
				right_subregion.box[cutfeat].first = cutval;
				right_subregion.depth = depth + 1;
				temp = divideByNameRatio(idx, right, right_subregion, depth + 1);
				node->child2 = temp;
				temp->parent = node;
				node->sub.pivot = cutval;
			}
			return node;
		}

		void computeMinMax(IndexType *ind, IndexType count, int element, ElementType &min_elem, ElementType &max_elem) {
			min_elem = getDataset(ind[0], element);
			max_elem = getDataset(ind[0], element);
			for (IndexType i = 1; i < count; ++i) {
				ElementType val = getDataset(ind[i], element);
				if (val < min_elem) min_elem = val;
				if (val > max_elem) max_elem = val;
			}
		}

		void divisionByPoint(IndexType *ind, IndexType count, IndexType &index, int &cutfeat, ElementType &cutval,
			const BoundingBox &bbox, size_t depth)
		{
			cutfeat = depth % m_dim;
			// for a balanced kd-tree, split in the median value
			std::vector<IndexType> cur_idx(count);
			for (IndexType i = 0; i < count; ++i) {
				cur_idx[i] = ind[i];
			}
			std::nth_element(cur_idx.begin(), cur_idx.begin() + cur_idx.size() / 2, cur_idx.end(),
				[this, &cutfeat](const IndexType a, const IndexType b) {
					return this->getDataset(a, cutfeat) < this->getDataset(b, cutfeat);
				});
			ElementType split_val = getDataset(cur_idx[cur_idx.size() / 2], cutfeat);

			ElementType min_elem, max_elem;
			computeMinMax(ind, count, cutfeat, min_elem, max_elem);

			if (split_val < min_elem) cutval = min_elem;
			else if (split_val > max_elem) cutval = max_elem;
			else cutval = split_val;

			IndexType lim1, lim2;
			planeSplit(ind, count, cutfeat, cutval, lim1, lim2);

			if (lim1 > count / 2) index = lim1;
			else if (lim2 < count / 2) index = lim2;
			else index = count / 2;
		}

		void divisionByRatio(IndexType *ind, IndexType count, IndexType &index, int &cutfeat, ElementType &cutval,
			const BoundingBox &bbox, size_t depth)
		{
			ofec::Real sum1 = 0.0;
			ofec::Real sum2 = 0.0;
			cutfeat = (m_relative_depth + depth) % m_dim;
			//cutfeat = rand() % m_dim;
			std::vector<IndexType> cur_idx(count);
			for (IndexType i = 0; i < count; ++i)
			{
				cur_idx[i] = ind[i];
			}

			index = cur_idx[cur_idx.size() / 2];

			for (auto &i : cur_idx)
				sum1 += m_ratio_data[i];

			for (auto j = 0; j < cur_idx.size() / 2; ++j)
			{
				sum2 += m_ratio_data[cur_idx[j]];
			}

			cutval = bbox.box[cutfeat].first + (bbox.box[cutfeat].second - bbox.box[cutfeat].first) * (sum2 / sum1);
		}

		void divisionByNameRatio(IndexType *ind, IndexType count, IndexType &index, int &cutfeat, ElementType &cutval,
			const BoundingBox &bbox, size_t depth)
		{
			ofec::Real sum1 = 0.0;
			ofec::Real sum2 = 0.0;
			cutfeat = (m_relative_depth + depth) % m_dim;
			//cutfeat = rand() % m_dim;
			std::vector<IndexType> cur_idx(count);
			for (IndexType i = 0; i < count; ++i) {
				cur_idx[i] = ind[i];
			}

			index = cur_idx[cur_idx.size() / 2];

			for (auto &i : cur_idx) {
				sum1 += m_name_ratio_data[i].second;
			}
			for (auto j = 0; j < cur_idx.size() / 2; ++j) {
				sum2 += m_name_ratio_data[cur_idx[j]].second;
			}

			cutval = bbox.box[cutfeat].first + (bbox.box[cutfeat].second - bbox.box[cutfeat].first) * (sum2 / sum1);
		}

		/**
		*  Subdivide the list of points by a plane perpendicular on axe corresponding
		*  to the 'cutfeat' dimension at 'cutval' position.
		*
		*  On return:
		*  dataset[ind[0..lim1-1]][cutfeat]<cutval
		*  dataset[ind[lim1..lim2-1]][cutfeat]==cutval
		*  dataset[ind[lim2..count]][cutfeat]>cutval
		*/
		void planeSplit(IndexType *ind, const IndexType count, int cutfeat, ElementType cutval, IndexType &lim1, IndexType &lim2)
		{
			/* Move vector indices for left subtree to front of list. */
			IndexType left = 0;
			IndexType right = count - 1;
			for (;;) {
				while (left <= right && getDataset(ind[left], cutfeat) < cutval) ++left;
				while (right && left <= right && getDataset(ind[right], cutfeat) >= cutval) --right;
				if (left > right || !right) break;  // "!right" was added to support unsigned Index types
				std::swap(ind[left], ind[right]);
				++left;
				--right;
			}
			/* If either list is empty, it means that all remaining features
			* are identical. Split in the middle to maintain a balanced tree.
			*/
			lim1 = left;
			right = count - 1;
			for (;;) {
				while (left <= right && getDataset(ind[left], cutfeat) <= cutval) ++left;
				while (right && left <= right && getDataset(ind[right], cutfeat) > cutval) --right;
				if (left > right || !right) break;  // "!right" was added to support unsigned Index types
				std::swap(ind[left], ind[right]);
				++left;
				--right;
			}
			lim2 = left;
		}

		NodePtr enqury(const std::vector<ElementType> &p, NodePtr node) const {
			if (node->child1 == nullptr && node->child2 == nullptr) {
				return node;
			}
			if (m_mode == 1) {
				if (p[node->sub.divfeat] < getDataset(node->sub.idx_sample, node->sub.divfeat)) {
					return enqury(p, node->child1);
				}
				else {
					return enqury(p, node->child2);
				}
			}
			else if (m_mode == 2 || m_mode == 3) {
				if (p[node->sub.divfeat] < node->sub.pivot) {
					return enqury(p, node->child1);
				}
				else {
					return enqury(p, node->child2);
				}
			}
			return nullptr;
		}

		NodePtr enqury(const ElementType *p, NodePtr node) const {
			if (node->child1 == nullptr && node->child2 == nullptr) {
				return node;
			}
			if (m_mode == 1) {
				if (p[node->sub.divfeat] < getDataset(node->sub.idx_sample, node->sub.divfeat)) {
					return enqury(p, node->child1);
				}
				else {
					return enqury(p, node->child2);
				}
			}
			else if (m_mode == 2 || m_mode == 3) {
				if (p[node->sub.divfeat] < node->sub.pivot) {
					return enqury(p, node->child1);
				}
				else {
					return enqury(p, node->child2);
				}
			}
			return nullptr;
		}

		void leafParent(IndexType idx_region, NodePtr node, NodePtr parent, NodePtr &result) {
			if (node->child1 == nullptr && node->child2 == nullptr) {
				if (node->lr.idx_region == idx_region) {
					if (node != m_root) result = parent;
					else {
						result = m_root;
					}
				}
				return;
			}
			if (node->child1 != nullptr && result == nullptr)  leafParent(idx_region, node->child1, node, result);
			if (node->child2 != nullptr && result == nullptr)  leafParent(idx_region, node->child2, node, result);
		}

		void boxRatio(BoundingBox &it, unsigned idx) {
			it.rat = 1;
			for (int i = 0; i < m_dim; ++i) {
				it.rat *= (it.box[i].second - it.box[i].first) / (m_init_region.box[i].second - m_init_region.box[i].first);
			}
			if (it.rat > m_lrat) {
				m_lrat = it.rat;
				m_lbox = idx;
			}
			if (it.rat < m_srat) {
				m_srat = it.rat;
				m_sbox = idx;
			}

			it.volume = 1;
			for (int i = 0; i < m_dim; ++i)
				it.volume *= (it.box[i].second - it.box[i].first);
		}

		void leafNode(IndexType idx_region, NodePtr node, NodePtr &leafnode) {
			if (node->child1 == nullptr && node->child2 == nullptr && node->lr.idx_region == idx_region) {
				leafnode = node;
				return;
			}
			if (node->child1 != nullptr)  leafNode(idx_region, node->child1, leafnode);
			if (node->child2 != nullptr)  leafNode(idx_region, node->child2, leafnode);
		}

		void getLeafRegions(NodePtr node, std::vector<size_t> &idx_regions) const {
			if (node->child1 == nullptr && node->child2 == nullptr) {
				idx_regions.push_back(node->lr.idx_region);
				return;
			}
			if (node->child1 != nullptr)
				getLeafRegions(node->child1, idx_regions);
			if (node->child2 != nullptr)
				getLeafRegions(node->child2, idx_regions);
		}

		void findUnionsAtDepth(size_t depth, NodePtr node, std::vector<std::vector<IndexType>> &region_unions) const {
			if (node->depth == depth) {
				std::vector<IndexType> idx_regions;
				getLeafRegions(node, idx_regions);
				region_unions.push_back(std::move(idx_regions));
				return;
			}
			if (node->child1 != nullptr)
				findUnionsAtDepth(depth, node->child1, region_unions);
			if (node->child2 != nullptr)
				findUnionsAtDepth(depth, node->child2, region_unions);
		}

		void neighborCheck(IndexType idx, NodePtr node, std::list<IndexType> &neighbors) const {
			if (node->child1 == nullptr && node->child2 == nullptr && idx != node->lr.idx_region) {
				neighbors.emplace_back(node->lr.idx_region);
				return;
			}
			if (node->child1 != nullptr) {
				if (!(m_subregions[idx]->box[node->sub.divfeat].second < node->sub.low || node->sub.pivot < m_subregions[idx]->box[node->sub.divfeat].first))
					neighborCheck(idx, node->child1, neighbors);
			}
			if (node->child2 != nullptr) {
				if (!(m_subregions[idx]->box[node->sub.divfeat].second < node->sub.pivot || node->sub.high < m_subregions[idx]->box[node->sub.divfeat].first))
					neighborCheck(idx, node->child2, neighbors);
			}
		}

		void mergeNode(NodePtr node, std::list<IndexType> &merged_regs) {
			if (node == nullptr) return;
			if (node->child1 != nullptr)
				mergeNode(node->child1, merged_regs);
			NodePtr parent = node->parent;
			NodePtr sibling = parent->child1 == node ? parent->child2 : parent->child1;
			if (sibling->child1 != nullptr)
				mergeNode(sibling->child1, merged_regs);
			m_subregions[node->lr.idx_region]->volume += m_subregions[sibling->lr.idx_region]->volume;
			m_subregions[node->lr.idx_region]->depth--;
			m_subregions[node->lr.idx_region]->box[parent->sub.divfeat].first = std::min(m_subregions[node->lr.idx_region]->box[parent->sub.divfeat].first, m_subregions[sibling->lr.idx_region]->box[parent->sub.divfeat].first);
			m_subregions[node->lr.idx_region]->box[parent->sub.divfeat].second = std::max(m_subregions[node->lr.idx_region]->box[parent->sub.divfeat].second, m_subregions[sibling->lr.idx_region]->box[parent->sub.divfeat].second);
			m_size--;
			parent->lr.idx_region = node->lr.idx_region;
			m_shelved_nodes.push_back(parent->child1);
			m_shelved_nodes.push_back(parent->child2);
			parent->child1 = nullptr;
			parent->child2 = nullptr;
			merged_regs.push_back(sibling->lr.idx_region);
			m_shelved_indices.insert(sibling->lr.idx_region);
		}

		/* Copy leaf nodes of p2 to p1 */
		void copyLeafNodes(NodePtr p1, NodePtr p2, std::map<IndexType, IndexType> &add_id) {
			if (p2->child1 == nullptr && p2->child2 == nullptr) {
				m_name_to_leaf_node[p1->name = p2->name] = p1;
				if (m_shelved_indices.empty()) {
					p1->lr.idx_region = m_subregions.size();
					m_subregions.emplace_back(new BoundingBox);
				}
				else {
					p1->lr.idx_region = *m_shelved_indices.begin();
					m_shelved_indices.erase(m_shelved_indices.begin());
				}
				m_subregions[p1->lr.idx_region]->depth = p1->depth;
				m_subregions[p1->lr.idx_region]->box = p2->lr.region->box;
				boxRatio(*m_subregions[p1->lr.idx_region], p1->lr.idx_region);
				p1->lr.region = m_subregions[p1->lr.idx_region].get();
				m_subregions[p1->lr.idx_region]->name = p1->name;
				p1->child1 = p1->child2 = nullptr;
				add_id[p2->lr.idx_region] = p1->lr.idx_region;
				m_size++;
			}
			else {
				p1->sub.divfeat = p2->sub.divfeat;
				p1->sub.pivot = p2->sub.pivot;
				p1->sub.low = p2->sub.low;
				p1->sub.high = p2->sub.high;

				m_name_to_leaf_node.erase(p1->name);

				if (!m_shelved_nodes.empty()) {
					p1->child1 = m_shelved_nodes.front();
					m_shelved_nodes.pop_front();
				}
				else
					p1->child1 = m_pool.allocate();
				p1->child1->depth = p1->depth + 1;
				p1->child1->parent = p1;
				copyLeafNodes(p1->child1, p2->child1, add_id);

				if (!m_shelved_nodes.empty()) {
					p1->child2 = m_shelved_nodes.front();
					m_shelved_nodes.pop_front();
				}
				else
					p1->child2 = m_pool.allocate();
				p1->child2->depth = p1->depth + 1;
				p1->child2->parent = p1;
				copyLeafNodes(p1->child2, p2->child2, add_id);
			}
		}
	};
	/** @} */ // end of grouping
} // end of NS


#endif /* NANOFLANN_HPP_ */
