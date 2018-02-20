#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
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
    bool m_store;
    ATOOLS::Blob_List * p_bl;
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
      m_store = reader.GetValue<int>("STORE_IN_HEPMC", 0);
      msg_Debugging()<<METHOD<<"(): { mode \n" <<
                       "skip b's from hard decay handler:  modus(" << m_modus <<
                       "), NLO(" << m_nlo << ") \n }" << std::endl;

      return true;
    }


    bool Run(ATOOLS::Blob_List *const bl){
      p_bl = bl;

      if (m_modus==0) return ReturnFunc(true);

      if (m_modus!=8){
          // go into signaif (NoB(ampl)) return true;l process blob and look for b-quarks
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
          ATOOLS::Blob *tmp(bl->FindLast(ATOOLS::btp::Shower));
          msg_Debugging() << *tmp << "\n";
          msg_Debugging() << "break" << std::endl;


          //try to acces shower amplitude
          // p_eventhandler: cant get p_phases!





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

                8: use cluster history for Veto

      */
          // msg_Info() << "Particles from DEC: " << num_dec << "    " << std::endl; //check

          if(m_modus==1 && ((numb_ps+numb_me)>0)) return ReturnFunc(false);
          if(m_modus==2 && numb_me>1) return ReturnFunc(false);
          if(m_modus==3 && numb_me>0) return ReturnFunc(false);
          if(m_modus==4){
              if(numb_me==1) return ReturnFunc(true);
              if (numb_me>0 || numb_ps>0) return ReturnFunc(false);
            }
          if(m_modus==5){
              if (numb_me>0 || numb_ps_all>0) return ReturnFunc(false);
            }
          if(m_modus==6){
              if (num_light_fs_quarks_me + num_fs_gluons_me >= 2) return ReturnFunc(true);
              if ((numb_ps+numb_me)>0) return ReturnFunc(false);

            }

          if(m_modus==7){
              /* still some ambiguities: There could be configurations in the ttjets event with a ME including two light jets and one b-quark.
           * If both light quarks are harder than the b-Quarks, this configuration should be kept in principle.
           * May one can check this in the cluster configuration?
           * One could for example directlx check the final cluster configuration and veto all configurations which have a ttbb cluster step?
           */
              if (numb_me>0) return ReturnFunc(false);
              if (m_nlo){
                  if (num_light_fs_quarks_me + num_fs_gluons_me > 2) return ReturnFunc(true);
                  if (num_light_fs_quarks_me >=2 ) return ReturnFunc(true);
                } else {
                  if (num_light_fs_quarks_me + num_fs_gluons_me > 1) return ReturnFunc(true);
                  if (num_light_fs_quarks_me >=1 ) return ReturnFunc(true);
                }
              if (numb_ps>0) return ReturnFunc(false);
            }

          return ReturnFunc(true);
        }

      else {
          /*   modus 8 -- use cluster history
                1) find first amplitude with b which does not come from a hard decay
                     -> no such amplitude: return true
                2) check, if previous amplitude has two light quarks or three lights partons in the FS
                    (check nlode mode for the exact number)
                     -> if yes: return true
                3) return false

         */

         ATOOLS::String_BlobDataBase_Map  bdmap = bl->FindFirst(ATOOLS::btp::Shower)->GetData();

         auto search = bdmap.find("AllAmplitudes");
         if (search==bdmap.end()) {
             THROW(fatal_error,"No matching amplitude found in blob. This algorithm works only with the ttbb224 branch!");
           }

         else{

         ATOOLS::Cluster_Amplitude * orig_ampl = search->second->Get<ATOOLS::Cluster_Amplitude*>();

         ATOOLS::Cluster_Amplitude * ampl = orig_ampl;


         // find first amplitude with b
         size_t pos=0;
         while (true){
             if (ampl==NULL) break;
             pos = FindB(ampl);
             if (pos==0) {
                 ampl = ampl->Next();
               }
             else break;
           }

         if (pos==0) return ReturnFunc(true);
         else return CheckPrev(ampl->Prev());

        }

      }
    }


    bool Finish(){}


    bool CheckPrev(ATOOLS::Cluster_Amplitude *ampl){
      //return only true, if a sufficient high number of light partons is present in this amplitude's final state
      if (ampl==NULL) return ReturnFunc(false);
      size_t num_q(0);
      size_t num_g(0);
      for (int i=2; i<ampl->Legs().size(); i++){
          ATOOLS::Cluster_Leg *leg = ampl->Legs().at(i);
          if (leg->Flav().IsGluon() && !leg->FromDec()) num_g++;
          if (abs(leg->Flav().Kfcode())<5 && !leg->FromDec()) num_q++;
      }

      if (m_nlo){
          if (num_q >=2 || (num_q + num_g)>2) return ReturnFunc(true);
          else return ReturnFunc(false);
        }
      else{
          if ((num_q + num_g)>=1) return ReturnFunc(true);
          else return ReturnFunc(false);
        }
    }

    size_t FindB(ATOOLS::Cluster_Amplitude *ampl){
      if (ampl==NULL) return 0;
      for (int i=2; i<ampl->Legs().size(); i++){
          ATOOLS::Cluster_Leg *leg = ampl->Legs().at(i);
          if ( abs(leg->Flav().Kfcode())==5 && !leg->FromDec()) return i;
        }
      return 0;
    }



    bool ReturnFunc(bool retval){
      if (retval==false){
          msg_Debugging() << "false\n";
        }

      if (!m_store) return retval;
      p_bl->FindFirst(ATOOLS::btp::Signal_Process)->AddData("Veto",new ATOOLS::Blob_Data<int>(!retval));
      return true;

    }



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

