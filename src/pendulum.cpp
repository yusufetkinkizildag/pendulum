#include <iostream>
#include <cmath>
#include <array>

namespace utility
{
    constexpr static auto G{9.80665}; // gravity
    constexpr static auto PI{3.14159265359};
    constexpr static auto T{0.000001}; // small time step
    constexpr static auto rad_to_deg{[](auto const r) noexcept{return r * 180 / PI;}};
} // namespace utility


class Pendulum
{
private:
    double theta; // current angle in radians
    double dtheta_dt; // first derivative of the theta wrt time
    double m; // mass
    double l; // length
    double k;
public:
    Pendulum(double const t_theta, double const t_dtheta_dt, double const t_m, double const t_l, double const t_k) noexcept
        : theta{t_theta}, dtheta_dt{t_dtheta_dt}, m{t_m}, l{t_l}, k{t_k}
    {
    }
    double operator()(double const tau, double const T) noexcept
    {
        auto const ddtheta_ddt{(tau - m * utility::G * l * std::sin(theta) - k * dtheta_dt) / (m * l * l)};
        dtheta_dt += ddtheta_ddt * T;
        theta += dtheta_dt * T;
        return utility::rad_to_deg(theta);
    };
};

int main()
{
    constexpr auto m{0.5};
    constexpr auto l{1.0};
    constexpr auto k{0.5};
    auto t{0.0};
    Pendulum pendulum{0.0, 0.0, m, l, k};
    std::array<double, 10> const taus{{0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0}};
    for (auto const tau : taus)
    {
        auto const theta_deg{pendulum(tau, t)};
        t += utility::T;
        std::cout << t << ',' << theta_deg << ',' << tau << '\n';
    }
    return 0;
}
