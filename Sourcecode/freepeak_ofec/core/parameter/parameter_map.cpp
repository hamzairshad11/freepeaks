#include "parameter_map.h"

namespace ofec {
    //bool operator==(const ParameterMap& p1, const ParameterMap& p2)
    //{
    //	return p1.m_map == p2.m_map;
    //}

    ParameterVariant& ParameterMap::operator[](const std::string& key) {
        return m_map[key];
    }

    const ParameterVariant& ParameterMap::at(const std::string& key) const {
        return m_map.at(key);
    }

    bool ParameterMap::has(const std::string& key) const {
        return m_map.count(key) > 0;
    }

    size_t ParameterMap::erase(const std::string& key) {
        return m_map.erase(key);
    }

    void ParameterMap::clear() {
        m_map.clear();
    }


    inline std::ostream& operator<<(std::ostream& os, const ParameterMap& pm)
    {
        os << "ParameterMap_begin\n";  // 输出开始标记


        for (const auto& [key, value] : pm.m_map) {
            os << key << " ";

            std::visit([&](const auto& v) {
                using T = std::decay_t<decltype(v)>;

                if constexpr (std::is_same_v<T, bool>) {
                    os << toString(ParameterType::kBool) << " " << v;
                }
                else if constexpr (std::is_same_v<T, int>) {
                    os << toString(ParameterType::kInt) << " " << v;
                }
                else if constexpr (std::is_same_v<T, size_t>) {
                    os << toString(ParameterType::kSizeT) << " " << v;
                }
                else if constexpr (std::is_same_v<T, char>) {
                    os << toString(ParameterType::kChar) << " " << v;
                }
                else if constexpr (std::is_same_v<T, Real>) {
                    os << toString(ParameterType::kReal) << " " << v;
                }
                else if constexpr (std::is_same_v<T, std::string>) {
                    os << toString(ParameterType::kString) << " " << v;
                }
                else if constexpr (std::is_same_v<T, std::vector<bool>>) {
                    os << toString(ParameterType::kVectorBool) << " " << v.size();
                    for (bool x : v) os << " " << x;
                }
                else if constexpr (std::is_same_v<T, std::vector<int>>) {
                    os << toString(ParameterType::kVectorInt) << " " << v.size();
                    for (int x : v) os << " " << x;
                }
                else if constexpr (std::is_same_v<T, std::vector<size_t>>) {
                    os << toString(ParameterType::kVectorSizeT) << " " << v.size();
                    for (size_t x : v) os << " " << x;
                }
                else if constexpr (std::is_same_v<T, std::vector<Real>>) {
                    os << toString(ParameterType::kVectorReal) << " " << v.size();
                    for (Real x : v) os << " " << x;
                }
                }, value);

            os << '\n';
        }

        os << "ParameterMap_end\n";  // 输出结束标记
        return os;
    }

    inline std::istream& operator>>(std::istream& in, ParameterMap& pm)
    {
        auto& map = pm.m_map;

        std::string key, type_str;

        // 读入 begin 标记
        std::string begin_tag;
        in >> begin_tag;
        if (begin_tag != "ParameterMap_begin")
            throw std::runtime_error("Missing ParameterMap_begin tag");

        while (true) {
            // 尝试读取 key，如果读到 end 标记则停止
            if (!(in >> key))
                break;
            if (key == "ParameterMap_end")
                break;

            in >> type_str;
            ParameterType type = fromString(type_str);

            switch (type) {
            case ParameterType::kBool: {
                bool v; in >> v;
                map[key] = v;
                break;
            }
            case ParameterType::kInt: {
                int v; in >> v;
                map[key] = v;
                break;
            }
            case ParameterType::kSizeT: {
                size_t v; in >> v;
                map[key] = v;
                break;
            }
            case ParameterType::kChar: {
                char v; in >> v;
                map[key] = v;
                break;
            }
            case ParameterType::kReal: {
                Real v; in >> v;
                map[key] = v;
                break;
            }
            case ParameterType::kString: {
                std::string v; in >> v;
                map[key] = v;
                break;
            }
            case ParameterType::kVectorBool: {
                size_t n; in >> n;
                std::vector<int> v(n);
                for (size_t i = 0; i < n; ++i) in >> v[i];
                std::vector<bool> bv(n);
                for (size_t i = 0; i < n; ++i) bv[i] = v[i];
                map[key] = bv;
                break;
            }
            case ParameterType::kVectorInt: {
                size_t n; in >> n;
                std::vector<int> v(n);
                for (auto& x : v) in >> x;
                map[key] = v;
                break;
            }
            case ParameterType::kVectorSizeT: {
                size_t n; in >> n;
                std::vector<size_t> v(n);
                for (auto& x : v) in >> x;
                map[key] = v;
                break;
            }
            case ParameterType::kVectorReal: {
                size_t n; in >> n;
                std::vector<Real> v(n);
                for (auto& x : v) in >> x;
                map[key] = v;
                break;
            }
            default:
                throw Exception("Unsupported ParameterType");
            }
        }

        return in;
    }

    inline bool operator==(const ParameterMap& lhs, const ParameterMap& rhs)
    {
        // 1️⃣ size 不同，直接 false
        if (lhs.m_map.size() != rhs.m_map.size()) {
            return false;
        }

        // 2️⃣ key + value 逐一比较
        for (const auto& [key, lhs_val] : lhs.m_map) {

            auto it = rhs.m_map.find(key);
            if (it == rhs.m_map.end()) {
                return false; // rhs 缺 key
            }

            const auto& rhs_val = it->second;

            // 3️⃣ variant index 不同，类型不一致
            if (lhs_val.index() != rhs_val.index()) {
                return false;
            }

            // 4️⃣ 内容比较
            bool equal = std::visit(
                [&](const auto& lv) -> bool {
                    using T = std::decay_t<decltype(lv)>;

                    // rhs 一定是同类型（上面 index 已检查）
                    const T& rv = std::get<T>(rhs_val);

                    if constexpr (std::is_same_v<T, std::vector<bool>>) {
                        if (lv.size() != rv.size()) return false;
                        for (size_t i = 0; i < lv.size(); ++i) {
                            if (static_cast<bool>(lv[i]) != static_cast<bool>(rv[i]))
                                return false;
                        }
                        return true;
                    }
                    else {
                        return lv == rv;
                    }
                },
                lhs_val
            );

            if (!equal) {
                return false;
            }
        }

        return true;
    }


}