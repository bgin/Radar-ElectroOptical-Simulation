
#ifndef __GMS_VERSION_H__
#define __GMS_VERSION_H__




namespace gms{

	namespace common{


	// Program version

		typedef struct VersionInfo {
			
		    unsigned int m_VersionMajor;

			unsigned int m_VersionMinor;

			unsigned int m_VersionMicro;

			unsigned int m_FullVersion;

			const char *  m_ProjectName;

			const char *  m_AppName;

			const char *  m_AppPurpose;

			const char *  m_AppType;

			const char *  m_AppArch;

			const char *   m_AppOS;

			const char *   m_AppCPU;

			const char *  m_AppCreateDate;

			const char *  m_AppLastBuild;

			const char *  m_AppCompilerVer;
			
		} VersionInfo_t;

		const VersionInfo_t gVersionInfo = {
			 1U,
			 0U,
			 0U,
			 1000U,
			 "Guided Missile Simulation",
			 "GuidedMissileSim",
			 "Air-to-Air Guided Missile Simulation and Modelling",
			 "Executable",
			 "x64",
			 "Linux",
			 "Intel Core i7 4770 HQ",
			 "Date: 2019-09-27 Time: 18:44 PM GMT+2",
			 "Date: To be specified",
			 "ICPC: 15.0.2.179 Build 20150121"
		};

   }
}


#endif /*__GMS_VERSION_H__*/
