#ifndef OFEC_FREE_PEAKS_FACTORY_H
#define OFEC_FREE_PEAKS_FACTORY_H

#include <map>
#include <set>
#include <functional>
#include <iterator> 
#include "../../../../core/problem/problem.h"

namespace ofec {
	namespace free_peaks {
		template<typename Base>
		struct FactoryFP {
			template<typename T>
			struct Register_ {
				Register_(const std::string& key) {
					ms_register_names.insert(key);
					map_.emplace(
						key,
						[]() {
							return T::create();

						}
					);
				}
			};

			static unsigned size() {
				return map_.size();
			}

			static Base* produce(const std::string& key) {
				auto it = map_.find(key);
				if (it == map_.end())
					THROW("the key " + key + " is not exist in FP register");
				auto base = it->second();
				base->setRegisterName(it->first);
				return base;
			}

			static Base* produceRandom(unsigned id) {
				auto it(map_.begin());
				advance(it, id);
				auto base = it->second();
				base->setRegisterName(it->first);
				return base;
			}

			// randnum in [0,1)
			static Base* produceRandom(double randnum) {
				unsigned id(map_.size() * randnum);
				auto it(map_.begin());
				advance(it, id);
				auto base = it->second();
				base->setRegisterName(it->first);
				return base;
			}

			//static const std::map<std::string, std::function<Base*()>>& get() {
			//	return map_;
			//}

			static const std::set<std::string>& registerNames() {
				return ms_register_names;
			}


			FactoryFP() = default;
		private:
			FactoryFP& operator=(const FactoryFP&) = delete;
			FactoryFP& operator=(FactoryFP&&) = delete;
			FactoryFP(const FactoryFP&) = delete;
			FactoryFP(FactoryFP&&) = delete;
			static std::set<std::string> ms_register_names;
			static std::map<std::string, std::function<Base* ()>>  map_;
		};

		template<typename Base>
		std::map<std::string, std::function<Base* ()>> FactoryFP<Base>::map_;

		template<typename Base>
		std::set<std::string> FactoryFP<Base>::ms_register_names;

#define REGISTER_FP(Base, Derived, key) FactoryFP<Base>::Register_<Derived> reg_##Derived(key);
	}
}

#endif // !OFEC_FREE_PEAKS_FACTORY_H