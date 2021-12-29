/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__it2des
#define _nrn_initial _nrn_initial__it2des
#define nrn_cur _nrn_cur__it2des
#define _nrn_current _nrn_current__it2des
#define nrn_jacob _nrn_jacob__it2des
#define nrn_state _nrn_state__it2des
#define _net_receive _net_receive__it2des 
#define castate castate__it2des 
#define evaluate_fct evaluate_fct__it2des 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gcabar _p[0]
#define shift _p[1]
#define qm _p[2]
#define qh _p[3]
#define m_inf _p[4]
#define tau_m _p[5]
#define h_inf _p[6]
#define tau_h _p[7]
#define m _p[8]
#define h _p[9]
#define cai _p[10]
#define cao _p[11]
#define Dm _p[12]
#define Dh _p[13]
#define ica _p[14]
#define carev _p[15]
#define phi_m _p[16]
#define phi_h _p[17]
#define v _p[18]
#define _g _p[19]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_evaluate_fct(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_it2des", _hoc_setdata,
 "evaluate_fct_it2des", _hoc_evaluate_fct,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gcabar_it2des", "mho/cm2",
 "shift_it2des", "mV",
 "tau_m_it2des", "ms",
 "tau_h_it2des", "ms",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"it2des",
 "gcabar_it2des",
 "shift_it2des",
 "qm_it2des",
 "qh_it2des",
 0,
 "m_inf_it2des",
 "tau_m_it2des",
 "h_inf_it2des",
 "tau_h_it2des",
 0,
 "m_it2des",
 "h_it2des",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 20, _prop);
 	/*initialize range parameters*/
 	gcabar = 0.003;
 	shift = 0;
 	qm = 2.5;
 	qh = 2.5;
 	_prop->param = _p;
 	_prop->param_size = 20;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _destexhe_IT2_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 20, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 it2des /Users/joao/Research/Models/NetPyNE/M1/mod/x86_64/destexhe_IT2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96485.3;
 static double R = 8.3145;
static int _reset;
static char *modelname = "Low threshold calcium current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int evaluate_fct(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int castate(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   evaluate_fct ( _threadargscomma_ v ) ;
   Dm = ( m_inf - m ) / tau_m ;
   Dh = ( h_inf - h ) / tau_h ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 evaluate_fct ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_m )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_h )) ;
  return 0;
}
 /*END CVODE*/
 static int castate (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   evaluate_fct ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_m)))*(- ( ( ( m_inf ) ) / tau_m ) / ( ( ( ( - 1.0 ) ) ) / tau_m ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_h)))*(- ( ( ( h_inf ) ) / tau_h ) / ( ( ( ( - 1.0 ) ) ) / tau_h ) - h) ;
   }
  return 0;
}
 
static int  evaluate_fct ( _threadargsprotocomma_ double _lv ) {
   m_inf = 1.0 / ( 1.0 + exp ( - ( _lv + shift + 50.0 ) / 7.4 ) ) ;
   h_inf = 1.0 / ( 1.0 + exp ( ( _lv + shift + 78.0 ) / 5.0 ) ) ;
   tau_m = ( 3.0 + 1.0 / ( exp ( ( _lv + shift + 25.0 ) / 10.0 ) + exp ( - ( _lv + shift + 100.0 ) / 15.0 ) ) ) / phi_m ;
   tau_h = ( 85.0 + 1.0 / ( exp ( ( _lv + shift + 46.0 ) / 4.0 ) + exp ( - ( _lv + shift + 405.0 ) / 50.0 ) ) ) / phi_h ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   phi_m = pow( qm , ( ( celsius - 24.0 ) / 10.0 ) ) ;
   phi_h = pow( qh , ( ( celsius - 24.0 ) / 10.0 ) ) ;
   evaluate_fct ( _threadargscomma_ v ) ;
   m = m_inf ;
   h = h_inf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   carev = ( 1e3 ) * ( R * ( celsius + 273.15 ) ) / ( 2.0 * FARADAY ) * log ( cao / cai ) ;
   ica = gcabar * m * m * h * ( v - carev ) ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  cai = _ion_cai;
  cao = _ion_cao;
 {   castate(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/joao/Research/Models/NetPyNE/M1/mod/destexhe_IT2.mod";
static const char* nmodl_file_text = 
  "TITLE Low threshold calcium current\n"
  ":\n"
  ":   Ca++ current responsible for low threshold spikes (LTS)\n"
  ":   RETICULAR THALAMUS\n"
  ":   Differential equations\n"
  ":\n"
  ":   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.\n"
  ":   The kinetics is described by standard equations (NOT GHK)\n"
  ":   using a m2h format, according to the voltage-clamp data\n"
  ":   (whole cell patch clamp) of Huguenard & Prince, J Neurosci.\n"
  ":   12: 3804-3817, 1992.\n"
  ":\n"
  ":    - Kinetics adapted to fit the T-channel of reticular neuron\n"
  ":    - Time constant tau_h refitted from experimental data\n"
  ":    - shift parameter for screening charge\n"
  ":\n"
  ":   Model described in detail in:   \n"
  ":     Destexhe, A., Contreras, D., Steriade, M., Sejnowski, T.J. and\n"
  ":     Huguenard, J.R.  In vivo, in vitro and computational analysis of\n"
  ":     dendritic calcium currents in thalamic reticular neurons.\n"
  ":     Journal of Neuroscience 16: 169-185, 1996.\n"
  ":   See also:\n"
  ":     http://www.cnl.salk.edu/~alain\n"
  ":     http://cns.fmed.ulaval.ca\n"
  ":\n"
  ":   Written by Alain Destexhe, Salk Institute, Sept 18, 1992\n"
  ":\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX it2des\n"
  "	USEION ca READ cai, cao WRITE ica\n"
  "	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift, qm, qh\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mV) =	(millivolt)\n"
  "	(mA) =	(milliamp)\n"
  "	(mM) =	(millimolar)\n"
  "\n"
  "	FARADAY = (faraday) (coulomb)\n"
  "	R = (k-mole) (joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	v		(mV)\n"
  "	celsius	= 36	(degC)\n"
  ":	eca	= 120	(mV)\n"
  "	gcabar	= .003	(mho/cm2)\n"
  "	shift	= 0 	(mV)\n"
  "	cai	= 2.4e-4 (mM)		: adjusted for eca=120 mV\n"
  "	cao	= 2	(mM)\n"
  "	qm	= 2.5\n"
  "	qh 	= 2.5\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m h\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ica	(mA/cm2)\n"
  "	carev	(mV)\n"
  "	m_inf\n"
  "	tau_m	(ms)\n"
  "	h_inf\n"
  "	tau_h	(ms)\n"
  "	phi_m\n"
  "	phi_h\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE castate METHOD cnexp\n"
  "	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)\n"
  "	ica = gcabar * m*m*h * (v-carev)\n"
  "}\n"
  "\n"
  "DERIVATIVE castate {\n"
  "	evaluate_fct(v)\n"
  "\n"
  "	m' = (m_inf - m) / tau_m\n"
  "	h' = (h_inf - h) / tau_h\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "INITIAL {\n"
  ":\n"
  ":   Activation functions and kinetics were obtained from\n"
  ":   Huguenard & Prince, and were at 23-25 deg.\n"
  ":   Transformation to 36 deg using Q10\n"
  ":\n"
  "	phi_m = qm ^ ((celsius-24)/10)\n"
  "	phi_h = qh ^ ((celsius-24)/10)\n"
  "\n"
  "	evaluate_fct(v)\n"
  "	m = m_inf\n"
  "	h = h_inf\n"
  "}\n"
  "\n"
  "PROCEDURE evaluate_fct(v(mV)) { \n"
  ":\n"
  ":   Time constants were obtained from J. Huguenard\n"
  ":\n"
  "\n"
  "	m_inf = 1.0 / ( 1 + exp(-(v+shift+50)/7.4) )\n"
  "	h_inf = 1.0 / ( 1 + exp((v+shift+78)/5.0) )\n"
  "\n"
  "	tau_m = ( 3 + 1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phi_m\n"
  "	tau_h = ( 85 + 1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ) ) / phi_h\n"
  "}\n"
  "UNITSON\n"
  ;
#endif
