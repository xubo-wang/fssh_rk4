//#define EIGEN_USE_MKL_ALL
//#define EIGEN_VECTORIZE_SSE4_2

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <complex>
//#pragma GCC system_headers
#include <D:/eigen/Eigen/Core>
#include <D:/eigen/Eigen/Dense>

#define DELTA 1.0
#define TOTAL_TIME 100.0
#define ITER 20

using namespace std;
using namespace Eigen;

class Output
//For the Traj class simulate() to return to a result.
//State 0 refers to lower_state, 1 refers to a higher state;
//trans refers to transmittion, reflect refers to refrection.
{
 public:
    Output(int current_state, double position)
	{
		state = current_state;
		trans = (position > 0);
		reflect = !trans;
	}
    Output __init()
    {
		state = 0;
		trans = 0;
		reflect = 0;
	}
	int state = 0;
	int trans = 0;
	int reflect = 0;
};
class Record
//low refers to lower state, high refers to higher state.
{
 public:
	int low_trans;
	int high_trans;
	int low_reflect;
	int high_reflect;
	Record __init()
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
class SimpleAvoidedCrossing
//Tully's first model in the paper.
//It's Hamiltonian can be learned from the V.
{
 public:
	SimpleAvoidedCrossing(double x)
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
//These parameters are assigned in the paper.
	const double a = 0.01;
	const double b = 1.6;
	const double c = 0.005;
	const double d = 1.0;
};

class ElectronicState
{
 public:
	ElectronicState(Matrix2d H, Matrix2d dH)
	{
		V = H;
		dV = dH;
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
	double ComputeForce(int state)
	{
		MatrixXcd out;
		out = MatrixXcd::Zero(1, 1);
		MatrixXcd state_vector;
		state_vector = basis_sets.col(state);
		out = state_vector.transpose() * (dV * state_vector);
		return -out(0, 0).real();
	}
	double ComputePotential(int state)
	{
		MatrixXcd out;
		out = Eigen::MatrixXcd::Zero(1, 1);
		MatrixXcd state_vector;
		state_vector = basis_sets.col(state);
		out = state_vector.transpose() * (V * state_vector);
		return -out(0, 0).real();
	}
	double ComputeCoupling(int bra, int ket)
	{
		MatrixXcd out;
		out = MatrixXcd::Zero(2, 2);
		double dE;
		if (bra != ket)
		{
			out = basis_sets.col(bra).transpose() * (dV * basis_sets.col(ket));
			dE = fabs(energy(bra, bra).real() - energy(ket, ket).real());
		}
		return out(0, 0).real() / dE;
	}
	Matrix2d ComputeNACMatrix(double velocity)
	{
		Matrix2d out;
		out << 0.0, 0.0, 0.0, 0.0;
		int i, j;
		for (i = 0; i<2; i++)
		{
			for (j = 0; j < i; j++)
			{
				double dij = ComputeCoupling(i, j);
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

class Traj
{
 public:
    Traj(double velo)
	{
		velocity = velo;
		pos = 0.0;
		rho << 1.0, 0.0, 0.0, 0.0;
		t = 0;
	}
	ElectronicState ComputeElectronics()
	{
		SimpleAvoidedCrossing sac = SimpleAvoidedCrossing(pos);
		return ElectronicState(sac.V, sac.dV);
	}
	ElectronicState ComputeElectronics(double pos)
	{
		SimpleAvoidedCrossing sac = SimpleAvoidedCrossing(pos);
		return ElectronicState(sac.V, sac.dV);
	}
	Output Simulate()
	{
		ElectronicState last_elec = ComputeElectronics();
		ElectronicState electronics = ComputeElectronics();
		
		double prob = 0.0;
		int step = 0;
		while (pos <= 5.0 && pos >= -5.0)
		{
			last_elec = electronics;
			ElectronicState es_1 = ComputeElectronics(pos);
			double a1 = electronics.ComputeForce(state) / mass;
			double v1 = dt * a1;

			double x1 = (v1 + velocity) * dt / 2 + pos;
			ElectronicState es_2 = ComputeElectronics(x1);
			double a2 = es_2.ComputeForce(state) / mass;
			double v2 = dt * a2;
			
			double x2 = (v2 + velocity) * dt / 2 + pos;
			ElectronicState es_3 = ComputeElectronics(x2);
			double a3 = es_3.ComputeForce(state) / mass;
			double v3 = dt * a3;

			double x3 = (v3 + velocity) * dt + pos;
			ElectronicState es_4 = ComputeElectronics(x3);
			double a4 = es_4.ComputeForce(state) / mass;
			double v4 = dt * a4;

			velocity += 1/3 * v1 + 1/6 * v2 + 1/6 * v3 + 1/3 * v4;
			pos += velocity * dt;
			electronics = ComputeElectronics();

			Propagate(electronics);
 			prob = SurfaceHopping(electronics);
			t += dt;
			energy = TotalEnergy(electronics);
			double potential = electronics.ComputePotential(state);
		}
		Output output(state, pos);
		return output;
	}
    void Propagate(ElectronicState es, double dt)
	{
		double velo = 0.5*(last_velo + velocity);//Euler
		Matrix2d D = es.ComputeNACMatrix(velo);

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
				G(col, row) = G(col, row) + i * D(col, row);
			}
		}
		Matrix2cd rho_prime;
		rho_prime << 0.0, 0.0, 0.0, 0.0;
		rho_prime = G * rho - rho * G;
		for (int col = 0; col < 2; col++)
		{
			for (int row = 0; row < 2; row++)
			{
				rho_prime(col, row) = rho_prime(col, row) * i;
			}
		}
		rho = rho + rho_prime * dt;
	}
	void Propagate(ElectronicState es)
	{
		double velo = 0.5 * (last_velo + velocity);//verlet
		Matrix2d D = es.ComputeNACMatrix(velo);

		int nstates = 2;

		Matrix2cd G;
		G.Zero(2, 2);
		G(0, 0) = es.energy(0, 0);
		G(1, 1) = es.energy(1, 1);
		complex<double> i(0, 1.0);
		for (int col = 0; col < 2; col++)
		{
			for (int row = 0; row < 2; row++)
			{
				G(col, row) = G(col, row) + i * D(col, row);
			}
		}
		Matrix2cd rho_prime;
		rho_prime << 0.0, 0.0, 0.0, 0.0;
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
	double SurfaceHopping(ElectronicState es)
	{
		int nstates = 2;

		double velo = 0.5 * (last_velo + velocity);
		Matrix2d W = es.ComputeNACMatrix(velo);

		vector<double> probs;

		int target = 0;
		for (; target < nstates; target++)
		{
			if (target == state)
				probs.push_back(0.0);
			else
			{
				double bij = -2.0 * rho(state, target).real() * W(target, state);
				if (bij < 0) bij = 0;
				probs.push_back(dt*bij / rho(state, state).real());
			}
		}
		double P = Sum(probs);//sum of probs
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
					double new_potential = es.ComputePotential(target);
					double old_potential = es.ComputePotential(state);
					double del_V = new_potential - old_potential;
					double coupling = es.ComputeCoupling(target, state);
					double component_kinetic = ModeKineticEnergy(coupling);
					if (del_V <= component_kinetic)
					{
						state = target;
						RescaleComponent(coupling, -del_V);
					}
				}
				break;
			}
		}
		return P;
	}
	double Sum(vector<double> vec) {
		double sum = 0;
		vector<double>::const_iterator it;
		for (it = vec.begin(); it != vec.end(); it++)
			sum += *it;
		return sum;
	}
	double random()
	{
		std::mt19937 gen(9527);
		std::uniform_real_distribution<> dis(0, 1);
		return dis(gen);
	}
    	double KineticEnergy()
	{
		return mass * pow(velocity, 2) / 2;
	}
	double TotalEnergy(ElectronicState es)
	{
		return KineticEnergy() + es.ComputePotential(state);
	}
	void RescaleComponent(double direction, double reduction)
	{
		direction /= fabs(direction);
		double Md = mass * direction;
		double a = Md * Md;
		double b = 2.0 * mass * velocity * Md;
		double c = -2.0 * mass * reduction;
		double scal = MinRoot(a, b, c);
		velocity += scal * direction;
	}
	double MinRoot(double a, double b, double c)
	{
		double x_1 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
		double x_2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
		return fabs(x_1) < fabs(x_2) ? x_1 : x_2;
	}

	double ModeKineticEnergy(double direction)
	{
		return 0.5 * mass * pow(velocity, 2);
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