#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <algorithm>

/*compile with:
  SHERPA_PREFIX=/home/s0118321/software/Sherpa/rel-2-2-4
  g++ -shared -g -I`$SHERPA_PREFIX/bin/Sherpa-config --incdir`  `$SHERPA_PREFIX/bin/Sherpa-config --ldflags`  -fPIC -o libSELANA.so SELANA.C
 */

namespace SELAN{


  class SELANA: public SHERPA::Analysis_Interface {

  protected:

  private:
    size_t m_modus;
    bool m_nlo;
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
      ATOOLS::Data_Reader reader(" ",";","//","=");
      reader.AddWordSeparator("\t");
      reader.SetAddCommandLine(false);
      reader.SetInputPath(m_inpath);
      std::string infile(m_infile);
      if (m_infile.find('|')!=std::string::npos)
        infile=infile.substr(0,infile.find('|'));
      reader.SetInputFile(infile+"|BEGIN_SELANA|END_SELANA");
      reader.SetComment("#");
      m_modus = reader.GetValue<int>("SELMODUS", 1);
      m_nlo = reader.GetValue<int>("NLO_Mode", 1);
      msg_Debugging()<<METHOD<<"(): { mode \n" <<
                       "skip b's from hard decay handler:  modus(" << m_modus <<
                       "), NLO(" << m_nlo << ") \n }" << std::endl;

      return true;
    }

    bool Run(ATOOLS::Blob_List *const bl){

      if (m_modus==0) return true;

      // go into signal process blob and look for b-quarks
      ATOOLS::Blob *sp(bl->FindFirst(ATOOLS::btp::Signal_Process));
      ATOOLS::Particle_Vector outvec(sp->GetOutParticles());
      size_t numb_me(0);
      size_t num_light_fs_quarks_me(0);
      size_t num_fs_gluons_me(0);

      for (size_t i=0; i<outvec.size(); i++){
          ATOOLS::Particle * particle(outvec.at(i));
          if (abs(particle->Flav().Kfcode())==5) {
              numb_me++;
            }
          if (abs(particle->Flav().Kfcode())<5) num_light_fs_quarks_me++;
          if (particle->Flav().IsGluon()) num_fs_gluons_me++;

        }
      // count number of outgoing b-quarks from the shower which do not go into the hard interaction
      ATOOLS::Blob *sh(bl->FindFirst(ATOOLS::btp::Shower));
      outvec=sh->GetOutParticles();
      size_t numb_ps(0), numb_ps_all(0);
      //size_t num_dec(0);
      for (size_t i(0); i<outvec.size();i++){
          ATOOLS::Particle * particle(outvec.at(i));
          //if(particle->FromDec()) num_dec++;
          if ( (abs(particle->Flav().Kfcode())==5) && !particle->FromDec() && particle->Info()!='G'){
              numb_ps++;
            }
          if ( abs(particle->Flav().Kfcode())==5 && particle->Info()!='G'){
              numb_ps_all++;
            }
        }

      /*
         modi:  0: no Veto, all events pass
                1: default, veto all b emissions eiter from ME or PS but not from HDH
                2: veto only ME-em with > 2 b's in ME
                3: veto all b ME-em
                4: veto Shower and ME, but not the "1 b from ME" case
                5: veto all all b emissions eiter from ME or PS, also from HDH
                6: veto all b emissions eiter from ME or PS but not from HDH,
                        exception: keep it, if >= 2 light jets in ME. (ttbb MC@NLO can only produce 1 extra light jet)
                7: veto all b emissions eiter from ME and all from PS but not from HDH,
                        exception: > 2 light jets from ME or >= 2 light quarks from ME

      */
     // msg_Info() << "Particles from DEC: " << num_dec << "    " << std::endl; //check

      if(m_modus==1 && ((numb_ps+numb_me)>0)) return false;
      if(m_modus==2 && numb_me>1) return false;
      if(m_modus==3 && numb_me>0) return false;
      if(m_modus==4){
          if(numb_me==1) return true;
          if (numb_me>0 || numb_ps>0) return false;
      }
      if(m_modus==5){
          if (numb_me>0 || numb_ps_all>0) return false;
      }
      if(m_modus==6){
          if (num_light_fs_quarks_me + num_fs_gluons_me >= 2) return true;
          if ((numb_ps+numb_me)>0) return false;

        }

      if(m_modus==7){
          /* still some ambiguities: There could be configurations in the ttjets event with a ME including two light jets and one b-quark.
           * If both light quarks are harder than the b-Quarks, this configuration should be kept in principle.
           * May one can check this in the cluster configuration?
           * One could for example directlx check the final cluster configuration and veto all configurations which have a ttbb cluster step?
           */
          if (numb_me>0) return false;
          if (m_nlo){
            if (num_light_fs_quarks_me + num_fs_gluons_me > 2) return true;
            if (num_light_fs_quarks_me >=2 ) return true;
          } else {
            if (num_light_fs_quarks_me + num_fs_gluons_me > 1) return true;
            if (num_light_fs_quarks_me >=1 ) return true;
          }
          if (numb_ps>0) return false;
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

