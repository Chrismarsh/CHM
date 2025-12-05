#include <concepts>

template<typename T, typename Data>
concept GuaranteeImplementExecute = requires(const T& t, Data& d) {
    { t.execute_impl(d) } -> std::same_as<void>;
};

template<class Derived, class Data>
class base_step {
protected:
    base_step() {};
public:
    void execute(Data& d) const {
        static_assert(GuaranteeImplementExecute<Derived,Data>,
                "Derived class must implement: void execute_impl(Data&) const");
        static_cast<const Derived*>(this)->execute_impl(d);
    }
    ~base_step() {};
};
