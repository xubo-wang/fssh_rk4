#include "fssh_rk4.cc"

int main()
{
	double k;
	vector<Record> P;
	Record single_k;
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
			double velo = (k);
			Traj fssh(velo);
			Output result = fssh.Simulate();
			single_k.low_reflect += (!result.state)*result.reflect;
			single_k.low_trans += (!result.state)*result.trans;
			single_k.high_reflect += result.reflect * result.state;
			single_k.high_trans += result.state*result.trans;
		}
		P.push_back(single_k);
	}
	vector<Record>::const_iterator it;
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


