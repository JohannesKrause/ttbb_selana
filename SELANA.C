#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Vector.H"
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
    size_t num_is_vetos, num_fs_vetos, num_both_vetos;

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
      m_modus = reader.GetValue<int>("SELMODUS", 8);
      m_nlo = reader.GetValue<int>("NLO_Mode", 1);
      m_store = reader.GetValue<int>("STORE_IN_HEPMC", 0);
      msg_Debugging()<<METHOD<<"(): { mode \n" <<
                       "skip b's from hard decay handler:  modus(" << m_modus <<
                       "), NLO(" << m_nlo << ") \n }" << std::endl;

      num_both_vetos=0;
      num_fs_vetos=0;
      num_is_vetos=0;

      return true;
    }


    bool Run(ATOOLS::Blob_List *const bl){
      p_bl = bl;

      if (m_modus==0) return ReturnFunc(true);

      if (m_modus!=8){
          //removed
          THROW(fatal_error,"Only modus 8 is implemented!");
        }


      else {
          /*   modus 8 -- use cluster history
           * final state configs:
                1) find first amplitude with b which does not come from a hard decay
                     -> no such amplitude: keep
                2) check, if previous amplitude has two light quarks or three lights partons in the FS
                    (check nlo mode for the exact number)
                     -> if yes: keep
                     -> otherwise:  do veto, contribution included in Xbb

            * initial state configs:
                1) check if previous amplitude has sufficient light final state partons -> keep
                2) count number of intermediate emissions, before initial-state b-splittings are unfolded
                  -> sufficient emissions? -> keep
                  -> if not: this is already part of the Xbb simulation and has to be removed

          do Veto, except both methods say "noVeto"

          It looks like the inital Veto is already included in the final Veto, all problematic configurations of IS-2) are found by the Final Veto, too.
          However, it is still usefull for the PDF-corrections
         */

         ATOOLS::String_BlobDataBase_Map  bdmap = bl->FindFirst(ATOOLS::btp::Shower)->GetData();

         auto search = bdmap.find("AllAmplitudes");
         if (search==bdmap.end()) {
             THROW(fatal_error,"No matching amplitude found in blob. This algorithm works only with the ttbb224 branch!");
           }

         else{

         ATOOLS::Cluster_Amplitude * orig_ampl = search->second->Get<ATOOLS::Cluster_Amplitude*>();

         ATOOLS::Cluster_Amplitude * ampl = orig_ampl;


         bool veto_final   = true;
         bool veto_initial = true;

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

         if (pos==0) veto_final=false;
         else veto_final = !CheckFinal(ampl->Prev());  // keep amplitude -> do not a veto

         //now do initial state treatment

         /* strategy for IS: look for first amplitude with b in IS
          *   no such amplitude:  allowed config
          * count number of IS-b's ?
          *   amplitude found:  check next and next->next if there was an IS splitting leading to final state b's
          *            -> check number of intermediate emissions. maximal one allowed (RS-emission)
          *         bb|x -> bb|xg -> bg| xgb  -> gg| xgbb     :  not allowed in tt+jets
          *         bb|x -> bb|xg -> bg| xgb  -> bg| xggb     :  allowed in tt+jets, two emissions harder than second b-splitting

         */


         // find first amplitude with b in IS and count number of b's there
         ampl = orig_ampl; // reset pointer
         size_t num_is_b=0;
         while (true){
             if (ampl==NULL) break;
             num_is_b = CountBinIS(ampl);
             if (num_is_b==0) {
                 ampl = ampl->Next();
               }
             else break;
           }
        //"ampl" points to amplitude with num_is_b b-quarks in IS
         if (num_is_b==0) veto_initial=false;
         else veto_initial = !CheckInitial(ampl);


         // fill statistics
         if (veto_final && veto_initial) num_both_vetos++;
         else{
             if (veto_final) num_fs_vetos++;
             if (veto_initial) num_is_vetos++;
           }


         /* ***** insert the counterterms for the PDF *******
          * This is only done, if veto_initial is true and if the b is in the initial state of an amplitude beeing maximally 2->3.
          * To all other amplitudes (2->4 or 2->5) the A2 terms would not contribute to the required order.
          *
          *  strategy:
          *    1) pick amplitude from above
          *
          *   starting with this amplitude, loop over all following:
          *
          *    2) check if amplitude has more than 3 FS-particles -> skip!
          *    3) for each new b in the initial state, calculate the correction weight.
          *
          * e.g.: a)  bb|tt -> bb|ttg: calculate correction for both b's for the first amplitude and
          *                            for one ratio, if there has been an initial state splitting?

          *       b)  gg|tt -> gb|ttb -> bb |ttbb:  calculate only correction for b-quark in second amplitude
          *       c)   gb|ttb -> xx|xxxx      : calculate only correction for first b-quark
          *
          * ******** questions:
          *
          *     * whats about addional IS-splittings on top off bb|tt?
          *       There the b-pdf enters the PDF-ratio in the Sudakovs, is this already higher order?
          *       -> Can one there just compute in the first amplitude the untouched b-pdf and in the secondstep the other one?
          *
          *  TODO:
          *      need Q2 in Amplitude. Which use? kT?
          *       was ist die richtige Skala? kTstart ist evtl höher als muf, falls es später einen Zerfall mit größerem kT gibt...
          *         -> nehme muF, falls 2->2 config, sonst kT_start für die Korrektur
          *
          *
          *   neue Strategie mache Korrektur für alle IS-b's bis zur 2->3 Konfig?
          */


         // start with the implementation

         // 2) if amplitude is to large, the corrections are higher order and not neccesary
         double weightfac(0.);
         int bookkeep[2] = {0,0};  // stores p_z of all b's which have already be corrected. the sign avoids double counting


         if(veto_initial){
             while(ampl->Legs().size()<=5){

                 for (int index=0; index < 2; index++){
                     if(abs(ampl->Leg(index)->Flav().Kfcode())==5) {
                         double test =  ampl->Leg(index)->Mom().operator [](3);
                         // if the ratio for one of the bookkeep's is >0, this has to be skipped
                         if ( !((bookkeep[0]/test >0) || (bookkeep[1]/test >0))  ) {
                             weightfac+=PdfCorrection(ampl,index);
                             bookkeep[index] =test;
                           }
                       }
                   }
                 ampl=ampl->Next();


               }
           }




         //do Veto, except both methods say "noVeto"
         return ReturnFunc(!(veto_final || veto_initial));

        }
      }
    }


    bool Finish(){

      //todo:  print some statistics: number of initial-vetos, number of final-vetos, combined vetoes?
      msg_Info() << "Some statistics about the veto procedure:" <<std::endl;
      msg_Info() << "   number of events with IS and FS veto: " << num_both_vetos << ".\n";
      msg_Info() << "   number of events with FS veto: " << num_fs_vetos << ".\n";
      msg_Info() << "   number of events with IS veto: " << num_is_vetos << ".\n";
    }

///////////////////////   helper functions final state configs   ///////////////////////////////////7

    bool CheckFinal(ATOOLS::Cluster_Amplitude *ampl){
      //return only true, if a sufficient high number of light partons is present in this amplitude's final state
      //  true:   keep amplitude
      //  false:  do veto
      if (ampl==NULL) return false; // could this be an allowed initial state configuration?
      size_t num_q(0);
      size_t num_g(0);
      for (int i=2; i<ampl->Legs().size(); i++){
          ATOOLS::Cluster_Leg *leg = ampl->Legs().at(i);
          if (leg->Flav().IsGluon() && !leg->FromDec()) num_g++;
          if (abs(leg->Flav().Kfcode())<5 && !leg->FromDec()) num_q++;
      }

      if (m_nlo){
          if (num_q >=2 || (num_q + num_g)>2) return true; // rethink
          else return false;
        }
      else{
          if ((num_q >=1) || (num_g>1)) return true; // rethink
          else return false;
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

    ///////////////////////   helper functions inital state configs   ///////////////////////////////////7


    size_t CountBinIS(ATOOLS::Cluster_Amplitude *ampl){
      if (ampl==NULL) return 0;
      size_t numb(0);
      for (int i=0; i<2; i++){
          ATOOLS::Cluster_Leg *leg = ampl->Legs().at(i);
          if ( abs(leg->Flav().Kfcode())==5 && !leg->FromDec()) numb++;
        }
      return numb;
    }

    size_t CountEW(ATOOLS::Cluster_Amplitude *ampl){
      if (ampl==NULL) return 0;
      size_t num(0);
      for (int i=2; i<ampl->Legs().size(); i++){
          ATOOLS::Cluster_Leg *leg = ampl->Legs().at(i);
          if ( !(leg->Flav().IsQuark() || leg->Flav().IsGluon()) ) num++; // todo: check, if this works.
        }
      return num;

    }

    bool CheckInitial(ATOOLS::Cluster_Amplitude *ampl){
      //return only true, if enough emissions take place before initial splitting is "unfolded"
      //  true:   keep amplitude
      //  false:  do veto
      if (ampl->Next()==NULL) return true;

      //      first: treatment of configs, where the initial b-splitting happens later in the evloution,
      //      after an already high enough number of emissions
      //    -> case, where the amplitude before the initial state splitting fullfills the number-of-light-jet condition

      if (ampl->Prev()) {
          if (CheckFinal(ampl->Prev())) return true;
        }

      //now: count number on intermediate emissions
      size_t num_intermed_emissions(0); // maximal 1 for nlo / 0 for LO
      size_t b_in_amp_prev=CountBinIS(ampl);
      size_t num_ew_particles_prev = CountEW(ampl);

      while (ampl->Next()){
          ampl= ampl->Next();
          size_t b_in_amp_new =CountBinIS(ampl);
          size_t num_ew_particles_new = CountEW(ampl);
          bool decay = (num_ew_particles_prev != num_ew_particles_new); // check, if this an actual emission or just a decay
          // identify decays, if number of non-QCD particles changes to the next amplitude
          if ((b_in_amp_new >= b_in_amp_prev) && !decay) num_intermed_emissions++;  // emission has not been unfolded in this step

          //check for num_intermed_emissions and decide if veto or not
          if (num_intermed_emissions > 1 && m_nlo) return true;
          if (num_intermed_emissions > 0 && !m_nlo) return true;

          // stop, if emission is unfolded
          if (b_in_amp_new==0) break;

          //reset counters
          b_in_amp_prev = b_in_amp_new;
          num_ew_particles_prev = num_ew_particles_new;

        }
      return false;

      }


    bool ReturnFunc(bool retval){
      if (retval==false){
          msg_Debugging() << "false\n";
          abort();
        }

      if (!m_store) return retval;
      p_bl->FindFirst(ATOOLS::btp::Signal_Process)->AddData("Veto",new ATOOLS::Blob_Data<int>(!retval));
      return true;

    }


    //  *****************************  PDF-Weight Correction ************************
    double PdfCorrection(ATOOLS::Cluster_Amplitude * ampl,int index){

      return 0.;
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

