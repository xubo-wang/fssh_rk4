//#define EIGEN_USE_MKL_ALL
//#define EIGEN_VECTORIZE_SSE4_2

#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <vector>
#include <D:/eigen/Eigen/Core>
#include <D:/eigen/Eigen/Dense>

#define DELTA 1.0
#define TOTAL_TIME 100.0
#define ITER 20

using namespace std;
using namespace Eigen;

class out_put
{
public:
	out_put(int _state, double position)
	{
		state = _state;
		trans = (position > 0);
		reflect = !trans;
	}
	out_put(int _i, int _j, int _k)
	{
		state = _i;
		trans = _j;
		reflect = _k;
	}
	int state = 0;
	int trans = 0;
	int reflect = 0;
};
class record
{
public:
	int low_trans;
	int high_trans;
	int low_reflect;
	int high_reflect;
	record __init()
	{
		low_trans = 0;
		high_trans = 0;
		low_reflect = 0;
		high_reflect = 0;
	}
	void clear()
	{
		low_trans = 0;
		high_trans = 0;
		low_reflect = 0;
		high_reflect = 0;
	}
};
class simple_avoided_crossing
{
public:
	simple_avoided_crossing(double x)
	{
		V(0, 0) = copysign(a, x) * (1.0 - exp(-b * fabs(x)));
		V(1, 1) = -1.0 * V(0, 0);
		V(0, 1) = c * exp(-d * x * x);
		V(1, 0) = V(0, 1);
		dV(0, 0) = a * b * exp(-b * fabs(x));
		dV(1, 1) = -dV(0, 0);
		dV(0, 1) = -2.0 * c * d * x * exp(-d * x * x);
		dV(1, 0) = dV(0, 1);
	}
	Matrix2d V;
	Matrix2d dV;
private:
	const double a = 0.01;
	const double b = 1.6;
	const double c = 0.005;
	const double d = 1.0;
};

class elec_state
{
public:
	elec_state(Matrix2d H, Matrix2d nabla_H)
	{
		V = H;
		dV = nabla_H;
		EigenSolver<Matrix2d> es(V);
		energy = es.eigenvalues().asDiagonal();
		basis_sets = es.eigenvectors();
	}
	int nstates()
	{
		return 2;
	}
	int ndim()
	{
		return 1;
	}
	double compute_force(int state)
	{
		MatrixXcd out;
		out = MatrixXcd::Zero(1, 1);
		MatrixXcd state_vec;
		state_vec = basis_sets.col(state);
		out = state_vec.transpose() * (dV * state_vec);
		return -out(0, 0).real();
	}
	double compute_potential(int state)
	{
		MatrixXcd out;
		out = MatrixXcd::Zero(1, 1);
		MatrixXcd state_vec;
		state_vec = basis_sets.col(state);
		out = state_vec.transpose() * (V * state_vec);
		return -out(0, 0).real();
	}
	double compute_coupling(int bra, int ket)
	{
		MatrixXcd out;
		out = MatrixXcd::Zero(1, 1);
		double dE;
		if (bra != ket)
		{
			out = basis_sets.col(bra).transpose() *(dV*basis_sets.col(ket));
			dE = fabs(energy(bra, bra).real() - energy(ket, ket).real());
		}
		return out(0, 0).real() / dE;
	}
	Matrix2d compute_NAC_matrix(double velocity)
	{
		Matrix2d out = Matrix2d::Zero(2, 2);
		int i, j;
		for (i = 0; i<2; i++)
		{
			for (j = 0; j < i; j++)
			{
				double dij = compute_coupling(i, j);
				out(i, j) = velocity * dij;
				out(j, i) = -out(i, j);
			}
		}
		return out;
	}
	Matrix2d V;
	Matrix2d dV;
	Matrix2cd energy;
	Matrix2cd basis_sets;
};

class traj
{
public:
	elec_state compute_elec()
	{
		simple_avoided_crossing sac = simple_avoided_crossing(pos);
		return elec_state(sac.V, sac.dV);
	}
	elec_state compute_elec(double pos)
	{
		simple_avoided_crossing sac = simple_avoided_crossing(pos);
		return elec_state(sac.V, sac.dV);
	}
	traj(double velo)
	{
		velocity = velo;
		pos = 0.0;
		rho << 1.0, 0.0, 0.0, 0.0;
		t = 0;
	}
	double kinetic_energy()
	{
		return mass * pow(velocity, 2) / 2;
	}
	double total_energy(elec_state es)
	{
		return kinetic_energy() + es.compute_potential(state);
	}
	out_put simulate()
	{
		elec_state last_elec = compute_elec();
		elec_state electronics = compute_elec();

		double ini_acc = electronics.compute_force(state) / mass;
		double velo = velocity;
		double dv = ini_acc * dt;
		double acc;

		double potential = electronics.compute_potential(state);
		double prob = 0.0;
		int step = 0;
		while (pos <= 5.0 && pos >= -5.0)
		{
			last_elec = electronics;
			electronics = compute_elec();
			double acc_1 = electronics.compute_force(state) / mass;
			double pos_1 = (acc_1 * dt + velocity) * dt;

			elec_state es_2 = compute_elec(pos + pos_1 / 2);
			double acc_2 = es_2.compute_force(state) / mass;
			double pos_2 = (acc_2*dt + velocity)*dt;

			elec_state es_3 = compute_elec(pos + pos_2 / 2);
			double acc_3 = es_3.compute_force(state) / mass;
			double pos_3 = (acc_3*dt + velocity)*dt;

			elec_state es_4 = compute_elec(pos + pos_3);
			double acc_4 = es_4.compute_force(state) / mass;

			velocity += 1 / 3 * acc_1*dt + 1 / 6 * acc_2 + 1 / 6 * acc_3 + 1 / 3 * acc_4;
			pos += velocity*dt;

			propagate(electronics);
			prob = surface_hopping(electronics, acc_1, acc_2, acc_3, acc_4);
			t += dt;
			energy = total_energy(electronics);
			potential = electronics.compute_potential(state);
		}
		out_put op = out_put::out_put(state, pos);
		return op;
	}
	double surface_hopping(elec_state es,double acc_1, double acc_2,double acc_3,double acc_4)//a little hesitate here about the rk4
	{
		int nstates = 2;

		double velo = 0.5 * (last_velo + velocity);
		Matrix2d W = es.compute_NAC_matrix(velo);

		vector<double> probs;

		int target = 0;
		for (; target < nstates; target++)
		{
			if (target == state)
				probs.push_back(0.0);
			else
			{
				double bij = -2.0 * rho(state, target).real() * W(target, state);
				probs.push_back(dt*bij / rho(state, state).real());
			}
		}
		double P = sum(probs);//sum of probs
		double accumulated_P = 0.0;
		double zeta = random();
		if (zeta < P)
		{
			for (target = 0; target < nstates; target++)
			{
				if (target == state)
				{
					continue;
				}
				accumulated_P += probs[target];
				if (zeta < accumulated_P)
				{
					double new_potential = es.compute_potential(target);
					double old_potential = es.compute_potential(state);
					double del_V = new_potential - old_potential;
					double coupling = es.compute_coupling(target, state);
					double component_kinetic = mode_kinetic_energy(coupling);
					if (del_V <= component_kinetic)
					{
						state = target;
						rescale_component(coupling, -del_V);
					}
				}
				break;
			}
		}
		return P;
	}
	double random()
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0, 1);
		return dis(gen);
	}
	double sum(vector<double> vec)
	{
		double sum = 0;
		vector<double>::const_iterator it;
		for (it = vec.begin(); it != vec.end(); it++)
			sum += *it;
		return sum;
	}
	void rescale_component(double direction, double reduction)
	{
		direction /= fabs(direction);
		double Md = mass * direction;
		double a = Md * Md;
		double b = 2.0 * mass * velocity * Md;
		double c = -2.0 * mass * reduction;
		double scal = min_root(a, b, c);
		velocity += scal * direction;
	}
	double min_root(double a, double b, double c)
	{
		double x_1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
		double x_2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
		return fabs(x_1) < fabs(x_2) ? x_1 : x_2;
	}
	double mode_kinetic_energy(double direction)
	{
		return 0.5 * mass * pow(velocity, 2);
	}
	void propagate(elec_state es, double dt)
	{
		double velo = 0.5*(last_velo + velocity);//Euler
		Matrix2d D = es.compute_NAC_matrix(velo);

		int nstates = 2;

		Matrix2cd G;
		G = Matrix2cd::Zero(2, 2);
		G(0, 0) = es.energy(0, 0);
		G(1, 1) = es.energy(1, 1);
		complex<double> i(0, 1.0);
		for (int col = 0; col < 2; col++)
		{
			for (int row = 0; row < 2; row++)
			{
				G(col, row) = G(col, row) + i*D(col, row);
			}
		}

		/*Matrix2cd diags;1
		Matrix2cd coeff;
		ComplexEigenSolver<Matrix2cd> ces(G);
		diags = ces.eigenvalues().asDiagonal();
		coeff = ces.eigenvectors();
		Matrix2cd coeff_t;
		coeff_t = coeff.conjugate();
		Matrix2cd tmp_rho = coeff_t * (rho * coeff);
		for (int i = 0; i < 2; i++)
		{
		for (int j = 0; j < 2; j++)
		{
		tmp_rho(i, j) = tmp_rho(i, j) * exp((diags(i) - diags(j))*dt);
		}
		}
		rho = coeff * (tmp_rho * coeff_t);*/
		Matrix2cd rho_prime = Matrix2cd::Zero(2, 2);
		rho_prime = G*rho - rho*G;
		for (int col = 0; col < 2; col++)
		{
			for (int row = 0; row < 2; row++)
			{
				rho_prime(col, row) = rho_prime(col, row) * i;
			}
		}

		rho = rho + rho_prime * dt;
	}
	void propagate(elec_state es)
	{
		double velo = 0.5*(last_velo + velocity);//Euler
		Matrix2d D = es.compute_NAC_matrix(velo);

		int nstates = 2;

		Matrix2cd G;
		G = Matrix2cd::Zero(2, 2);
		G(0, 0) = es.energy(0, 0);
		G(1, 1) = es.energy(1, 1);
		complex<double> i(0, 1.0);
		for (int col = 0; col < 2; col++)
		{
			for (int row = 0; row < 2; row++)
			{
				G(col, row) = G(col, row) + i*D(col, row);
			}
		}

		/*Matrix2cd diags;1
		Matrix2cd coeff;
		ComplexEigenSolver<Matrix2cd> ces(G);
		diags = ces.eigenvalues().asDiagonal();
		coeff = ces.eigenvectors();
		Matrix2cd coeff_t;
		coeff_t = coeff.conjugate();
		Matrix2cd tmp_rho = coeff_t * (rho * coeff);
		for (int i = 0; i < 2; i++)
		{
		for (int j = 0; j < 2; j++)
		{
		tmp_rho(i, j) = tmp_rho(i, j) * exp((diags(i) - diags(j))*dt);
		}
		}
		rho = coeff * (tmp_rho * coeff_t);*/
		Matrix2cd rho_prime = Matrix2cd::Zero(2, 2);
		rho_prime = G*rho - rho*G;
		for (int col = 0; col < 2; col++)
		{
			for (int row = 0; row < 2; row++)
			{
				rho_prime(col, row) = rho_prime(col, row) * i;
			}
		}

		rho = rho + rho_prime * dt;
	}
	const double mass = 2000.0;
	double velocity;
	double last_velo;
	double t;
	const double dt = DELTA;
	double pos = 0.0;
	int state = 0;
	const double nsteps = TOTAL_TIME / dt;
	Matrix2cd rho;
	double energy;
};
double gen_vel(double k)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<double> nd(k, k / 5 / sqrt(2));

	return nd(gen) / 2000;
}
int main()
{
	double k;
	vector<record> P;
	record single_k;
	ofstream file;
	file.open("fssh.log");
	double start = 1.0;
	double end = 40.0;
	double increment = 5.0;
	for (k = start; k < end; k += increment)
	{
		single_k.clear();
		for (int i = 0; i < ITER; i++)
		{
			double velo = gen_vel(k);
			traj fssh(velo);
			out_put result = fssh.simulate();
			single_k.low_reflect += (!result.state)*result.reflect;
			single_k.low_trans += (!result.state)*result.trans;
			single_k.high_reflect += result.reflect * result.state;
			single_k.high_trans += result.state*result.trans;
		}
		P.push_back(single_k);
	}
	vector<record>::const_iterator it;
	k = start;

	for (it = P.begin(); it != P.end(); it++)
	{
		int final_state = (*it).high_reflect + (*it).high_trans;
		int low_transmitt = (*it).low_trans;
		int low_reflect = (*it).low_reflect;
		int high_transmitt = (*it).high_trans;
		int high_reflection = (*it).high_reflect;
		file.setf(ios::right);
		file << k << ' ' << final_state << ' ' << low_transmitt << ' ' << low_reflect << ' ' << high_transmitt << ' ' << high_reflection << endl;
		k += increment;
	}

	file.close();
	return 0;
}

