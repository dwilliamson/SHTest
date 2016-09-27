
#include <sdla/RunApplication.h>
#include "SHTest.h"


int main(int argc, char* argv[])
{
	return (sdla::RunApplication<cSHTest>(1024, 768));
}