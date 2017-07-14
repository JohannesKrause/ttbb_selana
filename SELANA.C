#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <algorithm>

/*compile with:
  SHERPA_PREFIX=/home/s0118321/software/Sherpa/rel-2-2-2
  g++ -shared -g -I`$SHERPA_PREFIX/bin/Sherpa-config --incdir`  `$SHERPA_PREFIX/bin/Sherpa-config --ldflags`  -fPIC -o libSELANA.so SELANA.C
 */

namespace SELAN{


  class SELANA: public SHERPA::Analysis_Interface {

  protected:

  private:
    std::string m_inpath;
    std::string m_infile;
    std::string m_outpath;

  public:


    SELANA(const std::string &inpath,
           const std::string &infile,
           const std::string &outpath):SHERPA::Analysis_Interface("SELANA"),
      m_inpath(inpath), m_infile(infile), m_outpath(outpath)
    {
      msg_Debugging()<<"SELANA Selector activ."<<std::endl;
    }

    void ShowSyntax(const int i)
    {
      if (!msg_LevelIsInfo() || i==0) return;
      msg_Out()<<METHOD<<"(): {\n\n"
              <<"custom analysis acting as selector"
             <<"}"<<std::endl;
    }


    bool Init()
    {
      return true;
    }

    bool Run(ATOOLS::Blob_List *const bl){


      // veto b-emissions from ME
      ATOOLS::Blob *sp(bl->FindFirst(ATOOLS::btp::Signal_Process));
      ATOOLS::Particle_Vector outvec(sp->GetOutParticles());
      for (size_t i=0; i<outvec.size(); i++){
          ATOOLS::Particle * particle(outvec.at(i));
          if (abs(particle->Flav().Kfcode())==5) return false;
        }

      // veto b-emissions from the Shower which do not originate from the Hard Decay Handler
      ATOOLS::Blob *sh(bl->FindFirst(ATOOLS::btp::Shower));
      ATOOLS::Particle_Vector outvec(sh->GetOutParticles());
      for (size_t i(0); i<outvec.size();i++){
          ATOOLS::Particle * particle(outvec.at(i));
          if ( (abs(particle->Flav().Kfcode())==5) && !particle->Dec() ) return false;
        }

      return true;
    }

    bool Finish(){}

  };// end of class SELANA



} // end of namespace SELAN



using namespace SHERPA;
using namespace SELAN;

DECLARE_GETTER(SELANA,"SELANA",
               Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,SELANA>::
operator()(const Analysis_Arguments &args) const
{
  msg_Info()<<METHOD<<"(): {}"<<std::endl;
  return new SELANA(args.m_inpath,args.m_infile,args.m_outpath);
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,
SELANA>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Selector Analysis";
}

