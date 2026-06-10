/********* Begin Register Information **********
{
	"name": "free_peaks",
	"identifier": "FreePeaks",
	"tags": [ "continuous"]
}
*********** End Register Information **********/

#ifndef OFEC_FREE_PEAKS_H
#define OFEC_FREE_PEAKS_H

#include "../../../../core/problem/continuous/continuous.h"
#include "subspace_kdtree.h"
#include <mutex>

namespace ofec {
#define CAST_FPs(pro) dynamic_cast<FreePeaks*>(pro)

	class FreePeaks : virtual public Continuous {
		OFEC_CONCRETE_INSTANCE(FreePeaks)

	public:
		static const std::string ms_file_dir;
		static std::string directory();
		static void registerFP();


	protected:
		void addInputParameters();
		void initialize_(Environment* env) override;
		void evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons) const override;
		void calculateOptima();
		void updateOptima(Environment* env) override;
		void readFromFile(const std::string& file_path);
		void writeToFile() const;


		std::string m_generation_type, m_file_name;
		free_peaks::SubspaceKDTree m_subspace_tree;
		Real m_objective_accuracy;
		static bool ms_registered;
		static std::mutex ms_register_mtx;

	public:

		void updateCandidates(const SolutionBase& sol, std::vector<std::unique_ptr<SolutionBase>>& candidates) const;
		bool isSolved(const std::vector<std::unique_ptr<SolutionBase>>& candidates) const;
		std::vector<Real> errorToPeaks(const std::vector<std::unique_ptr<SolutionBase>>& candidates) const;


		const free_peaks::SubspaceKDTree& subspaceTree() const { return m_subspace_tree; }
		virtual free_peaks::Subproblem* subproblem(const std::string& subspace_name) const;
		virtual const std::vector<std::pair<Real, Real>>& subspaceBox(const std::string& subspace_name) const;
		std::string generationType() const { return m_generation_type; }
		std::string fileName() const { return m_file_name; }

		// free_peak  cosntrution functions
		void setSizes(size_t num_vars, size_t num_objs, size_t num_cons);
		// subtree_name  subspace_name , subspace ratio
		void setKDtree(const std::vector<std::pair<std::string, std::vector<std::pair<std::string, double>>>>& treeData);
		void setSubproblem(const std::string& subspace_name, free_peaks::Subproblem* subpro);
		void setSubproblem(const std::string& subspace_name, const ParameterMap& param);
		// bind data 
		void bindData();
		void outputTotalFile();

		virtual void recordInputParameters();
		virtual void restoreInputParameters();

	};
}

#endif // !OFEC_FREE_PEAKS_H