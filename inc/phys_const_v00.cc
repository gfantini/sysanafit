#include "phys_const_v00.hh"
TRandom3 rndmgen_pc;
PhysConst::PhysConst(std::string iname,Double_t ivalue,Double_t ierror,std::string iunits,std::string ireference="")
{
	name = iname;
	value = ivalue;
	error = ierror;
	units = iunits;
	reference = ireference;
}
void PhysConst::SetReference(std::string ireference)
{
	reference = ireference;
}

Double_t PhysConst::GetValue()
{
	return value;
}
Double_t PhysConst::GetError()
{
	return error;
}
std::string PhysConst::GetName()
{
	return name;
}
Double_t PhysConst::GetRandom()
{
	return rndmgen_pc.Gaus(value,error);
}
void PhysConst::PrintInfo()
{
	cerr << "[PhysConst::PrintInfo] value = " << value << endl;
	cerr << "[PhysConst::PrintInfo] error = " << error << endl;
	cerr << "[PhysConst::PrintInfo] units = " << units << endl;
	cerr << "[PhysConst::PrintInfo] reference = " << reference << endl;
}
