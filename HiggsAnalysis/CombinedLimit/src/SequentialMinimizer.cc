#include "../interface/SequentialMinimizer.h"

#include <cmath>
#include <memory>
#include <algorithm>
#include <limits>
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooLinkedListIter.h"
#include <Math/MinimizerOptions.h>
#include <boost/foreach.hpp>
#include "../interface/ProfilingTools.h"
#define foreach BOOST_FOREACH

//#define DEBUG_SMD
#ifdef DEBUG_SMD
#define printf if (RT_DEBUG_SMD) printf
#else
#define printf if (0) printf
#endif
namespace { 
    const double GOLD_R1 = 0.61803399 ;
    const double GOLD_R2 = 1-0.61803399 ;
    const double XTOL    = 10*std::sqrt(std::numeric_limits<double>::epsilon());
#ifdef DEBUG_SMD
    bool RT_DEBUG_SMD = 0;
#endif
}

bool OneDimMinimizer::minimize(int steps, double ytol, double xtol) 
{
    // get initial bracket
    seek();
    
    bool done = doloop(steps,ytol,xtol);
    parabolaStep();
    var_->setVal(xi_[1]);
    var_->setError(xi_[2] - xi_[0]);
    return done;
}

OneDimMinimizer::ImproveRet OneDimMinimizer::improve(int steps, double ytol, double xtol, bool force) 
{
    double x0 = var_->getVal();
    if (x0 < xi_[0] || x0 > xi_[2]) {
        // could happen if somebody outside this loop modifies some parameters
        printf("ODM: ALERT: variable %s outside bounds x = [%.4f, %.4f, %.4f], x0 = %.4f\n", var_->GetName(), xi_[0], xi_[1], xi_[2], x0);
        x0 = xi_[1]; 
        var_->setVal(x0);
    } else {
        xi_[1] = x0;
    }
    double y0 = nll_->getVal();
    yi_[1] = y0;
    yi_[0] = eval(xi_[0]);
    yi_[2] = eval(xi_[2]);
    if (xtol == 0) xtol = (xi_[1]+XTOL)*XTOL;

    printf("ODM: start of improve %s x = [%.4f, %.4f, %.4f], y = [%.4f, %.4f, %.4f]\n", var_->GetName(), xi_[0], xi_[1], xi_[2], yi_[0], yi_[1], yi_[2]);

    if (yi_[1] < yi_[0] && yi_[1] < yi_[2]) {
        if (ytol > 0 && (max(yi_[2],yi_[0]) - yi_[1]) < ytol) {
            printf("ODM: immediate ytol for %s: ymin %.8f, ymax %.8f, diff %.8f\n", var_->GetName(), yi_[1], max(yi_[2],yi_[0]), max(yi_[2],yi_[0]) - yi_[1]);
            if (!force || parabolaStep()) return Unchanged;
        }
        if (xtol > 0 && (xi_[2] - xi_[0]) < xtol) {
            printf("ODM: immediate xtol for %s: xmin %.8f, xmax %.8f, diff %.8f\n", var_->GetName(), xi_[0], xi_[2], xi_[2] - xi_[0]);
            if (!force || parabolaStep()) return Unchanged;
        }
    } else {
        reseek();
    }

    //post-condition: always a sorted interval
    assert(xi_[0] < xi_[2]);
    assert(xi_[0] <= xi_[1] && xi_[1] <= xi_[2]);
    // if midpoint is not not one of the extremes, it's not higher than that extreme
    assert(xi_[1] == xi_[0] || yi_[1] <= yi_[0]); 
    assert(xi_[1] == xi_[2] || yi_[1] <= yi_[2]);

    bool done = doloop(steps,ytol,xtol);
    parabolaStep(); 

    //post-condition: always a sorted interval
    assert(xi_[0] < xi_[2]);
    assert(xi_[0] <= xi_[1] && xi_[1] <= xi_[2]);
    // if midpoint is not not one of the extremes, it's not higher than that extreme
    assert(xi_[1] == xi_[0] || yi_[1] <= yi_[0]); 
    assert(xi_[1] == xi_[2] || yi_[1] <= yi_[2]);


    if (ytol > 0 && fabs(yi_[1] - y0) < ytol) {
        printf("ODM: final ytol for %s: ymin(old) %.8f, ymin(new) %.8f, diff %.8f\n", var_->GetName(), y0, yi_[1], y0 -yi_[1]);
        if (!force) var_->setVal(x0);
        return Unchanged;
    } 

    if (xtol > 0 && fabs(xi_[1] - x0) < xtol) {
        var_->setVal(force ? xi_[1] : x0);
        return Unchanged;
    }
    printf("ODM: doloop for %s is %s\n", var_->GetName(), done ? "Done" : "NotDone");
    printf("ODM: end of improve %s x = [%.4f, %.4f, %.4f], y = [%.4f, %.4f, %.4f]\n", var_->GetName(), xi_[0], xi_[1], xi_[2], yi_[0], yi_[1], yi_[2]);
    var_->setVal(xi_[1]);
    var_->setError(xi_[2] - xi_[0]);
    return done ? Done : NotDone; 
}

bool OneDimMinimizer::doloop(int steps, double ytol, double xtol) 
{
    if (steps <= 0) steps = 100;
    for (int i = 0; i < steps; ++i) {
        goldenBisection();
        if (ytol > 0 && (max(yi_[2],yi_[0]) - yi_[1]) < ytol) {
            printf("ODM: intermediate ytol for %s: ymin %.8f, ymax %.8f, diff %.8f\n", var_->GetName(), yi_[1], max(yi_[2],yi_[0]), max(yi_[2],yi_[0]) - yi_[1]);
            return true;
        }
        if (xtol > 0 && (xi_[2] - xi_[0]) < xtol) {
            return true;
        }
        printf("ODM: step %d/%d done for %s: ymin %.8f, ymax %.8f, diff %.8f\n", i+1, steps, var_->GetName(), yi_[1], max(yi_[2],yi_[0]), max(yi_[2],yi_[0]) - yi_[1]);
        printf("ODM: %s x = [%.4f, %.4f, %.4f], y = [%.4f, %.4f, %.4f]\n", var_->GetName(), xi_[0], xi_[1], xi_[2], yi_[0], yi_[1], yi_[2]);
    }

    return false;
}

void OneDimMinimizer::seek() 
{
    bool   hasmax = var_->hasMax();
    double xmin = var_->getMin(), xmax = var_->getMax();
    double xerr = hasmax ? 0.2*(xmax-xmin) : 1.0;
    if (var_->hasError() && var_->getError() > xerr) {
        xerr = var_->getError();
    }
    var_->setError(xerr);
    reseek();
}
void OneDimMinimizer::reseek() 
{
    double xmin = var_->getMin(), xmax = var_->getMax();
    double xerr = var_->getError(); 
    xi_[1] = var_->getVal(); 
    yi_[1] = eval(xi_[1]);
    xi_[0] = std::max(xmin, xi_[1]-xerr);
    yi_[0] = (xi_[0] == xi_[1] ? yi_[1] : eval(xi_[0]));
    xi_[2] = std::min(xmax, xi_[1]+xerr);
    yi_[2] = (xi_[2] == xi_[1] ? yi_[1] : eval(xi_[2]));
    assert(xi_[0] < xi_[2]);
    assert(xi_[0] <= xi_[1] && xi_[1] <= xi_[2]);

    for (;;) {
        //printf("ODM: bracketing %s in x = [%.4f, %.4f, %.4f], y = [%.4f, %.4f, %.4f]\n", var_->GetName(), xi_[0], xi_[1], xi_[2], yi_[0], yi_[1], yi_[2]);
        if (yi_[0] < yi_[1]) {
            assign(2,1); // 2:=1
            assign(1,0); // 1:=0
            xi_[0] = std::max(xmin, xi_[1]-xerr);
            yi_[0] = (xi_[0] == xi_[1] ? yi_[1] : eval(xi_[0]));
        } else if (yi_[2]  < yi_[1]) {
            assign(0,1); // 0:=1
            assign(1,2); // 1:=2
            xi_[2] = std::min(xmax, xi_[1]+xerr);
            yi_[2] = (xi_[2] == xi_[1] ? yi_[1] : eval(xi_[2]));
        } else {
            xerr /= 2;
            break;
        }
        xerr *= 2;
    }
    //printf("ODM: bracketed minimum of %s in [%.4f, %.4f, %.4f]\n", var_->GetName(), xi_[0], xi_[1], xi_[2]);
    //post-condition: always a sorted interval
    assert(xi_[0] < xi_[2]);
    assert(xi_[0] <= xi_[1] && xi_[1] <= xi_[2]);
    // if midpoint is not not one of the extremes, it's not higher than that extreme
    assert(xi_[1] == xi_[0] || yi_[1] <= yi_[0]); 
    assert(xi_[1] == xi_[2] || yi_[1] <= yi_[2]);
    var_->setError(xerr);
}

void OneDimMinimizer::goldenBisection() 
{
    //pre-condition: always a sorted interval
    assert(xi_[0] < xi_[2]);
    assert(xi_[0] <= xi_[1] && xi_[1] <= xi_[2]);
    // if midpoint is not not one of the extremes, it's not higher than that extreme
    assert(xi_[1] == xi_[0] || yi_[1] <= yi_[0]); 
    assert(xi_[1] == xi_[2] || yi_[1] <= yi_[2]);
    if (xi_[0] == xi_[1] || xi_[1] == xi_[2]) {
        int isame = (xi_[0] == xi_[1] ? 0 : 2);
        /// pre-condition: the endpoint equal to x1 is not the highest
        assert(yi_[isame] <= yi_[2-isame]);
        xi_[1] = 0.5*(xi_[0]+xi_[2]);
        yi_[1] = eval(xi_[1]);
        if (yi_[1] < yi_[isame]) {
            // maximum is in the interval-
            // leave as is, next bisection steps will find it
        } else {
            // maximum remains on the boundary, leave both points there
            assign(2-isame, 1);
            assign(1, isame); 
        }
    } else {
        int inear = 2, ifar = 0;
        if (xi_[2]-xi_[1] > xi_[1] - xi_[0]) {
            inear = 0, ifar = 2;
        } else {
            inear = 2, ifar = 0;
        }
        double xc = xi_[1]*GOLD_R1 + xi_[ifar]*GOLD_R2;
        double yc = eval(xc);
        //printf("ODM: goldenBisection:\n\t\tfar = (%.2f,%.8f)\n\t\tnear = (%.2f,%.8f)\n\t\tcenter  = (%.2f,%.8f)\n\t\tnew  = (%.2f,%.8f)\n",
        //            xi_[ifar],  yi_[ifar], xi_[inear], yi_[inear], xi_[1], yi_[1], xc, yc);
        if (yc < yi_[1]) {   // then use 1, c, far
            assign(inear, 1);
            xi_[1] = xc; yi_[1] = yc;
        } else {  // then use c, 1, near
            xi_[ifar] = xc; yi_[ifar] = yc;
        }
    }
    //post-condition: always a sorted interval
    assert(xi_[0] < xi_[2]);
    assert(xi_[0] <= xi_[1] && xi_[1] <= xi_[2]);
    // if midpoint is not not one of the extremes, it's not higher than that extreme
    assert(xi_[1] == xi_[0] || yi_[1] <= yi_[0]); 
    assert(xi_[1] == xi_[2] || yi_[1] <= yi_[2]);

}

double OneDimMinimizer::parabolaFit() 
{
    if (xi_[0] == xi_[1] || xi_[1] == xi_[2]) { 
        return xi_[1]; 
    }
    double dx0 = xi_[1] - xi_[0], dx2 = xi_[1] - xi_[2];
    double dy0 = yi_[1] - yi_[0], dy2 = yi_[1] - yi_[2];
    double num = dx0*dx0*dy2 - dx2*dx2*dy0;
    double den = dx0*dy2     - dx2*dy0;
    if (den != 0) {
        double x = xi_[1] - 0.5*num/den;
        if (xi_[0] < x && x < xi_[2]) {
            return x;
        }
    } 
    return xi_[1];
}

bool OneDimMinimizer::parabolaStep() {
    double xc = parabolaFit();
    if (xc != xi_[1]) {
        double yc = eval(xc);
        if (yc < yi_[1]) {
            xi_[1] = xc; 
            yi_[1] = yc;
            //post-condition: always a sorted interval
            assert(xi_[0] < xi_[2]);
            assert(xi_[0] <= xi_[1] && xi_[1] <= xi_[2]);
            // if midpoint is not not one of the extremes, it's not higher than that extreme
            assert(xi_[1] == xi_[0] || yi_[1] <= yi_[0]); 
            assert(xi_[1] == xi_[2] || yi_[1] <= yi_[2]);
            return true;
        }
    }
    return false;
}

double OneDimMinimizer::eval(double x) 
{
    double x0 = var_->getVal();
    var_->setVal(x); 
    double y = nll_->getVal();
    var_->setVal(x0);
    return y;
}

SequentialMinimizer::SequentialMinimizer(RooAbsReal *nll, RooRealVar *poi) :
    nll_(nll),
    state_(Cleared)
{
#ifdef DEBUG_SMD
    RT_DEBUG_SMD = runtimedef::get("DEBUG_SMD");
#endif
    std::auto_ptr<RooArgSet> args(nll->getParameters((const RooArgSet *)0));
    workers_.reserve(args->getSize());
    if (poi != 0) workers_.push_back(Worker(nll,poi));
    std::auto_ptr<TIterator> iter(args->createIterator());
    for (RooAbsArg *a = (RooAbsArg*) iter->Next(); a != 0; a = (RooAbsArg*) iter->Next()) {
        RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
        if (rrv != 0 && !rrv->isConstant() && rrv != poi) {
            workers_.push_back(Worker(nll,rrv));
        }
    }
}

bool SequentialMinimizer::minimize(double ytol, int bigsteps, int smallsteps) 
{
    foreach(Worker &w, workers_) {
        if (w.var().isConstant()) {
            w.state = Fixed;
        } else {
            w.minimize(1); 
            w.state = Ready; 
        }
    }
    state_ = Ready;
    return improve(ytol, bigsteps, smallsteps);
}

bool SequentialMinimizer::improve(double ytol, int bigsteps, int smallsteps)
{
    // catch improve before minimize case
    if (state_ == Cleared) return minimize(ytol,bigsteps,smallsteps);

    // setup default tolerances and steps
    if (ytol == 0) ytol = ROOT::Math::MinimizerOptions::DefaultTolerance()/sqrt(workers_.size());
    if (bigsteps == 0) bigsteps = 100 * (workers_.size());

    // list of done workers (latest-done on top)
    std::list<Worker*> doneWorkers;

    // start with active workers, for all except constants
    foreach(Worker &w, workers_) {
        if (w.var().isConstant()) w.state = Fixed;
        else w.state = Active;
    }

    state_ = Active;
    for (int i = 0; i < bigsteps; ++i) {
        printf("Start of loop. State is %s\n",(state_ == Done ? "DONE" : "ACTIVE"));
        State newstate = Done;
        foreach(Worker &w, workers_) {
            OneDimMinimizer::ImproveRet iret = OneDimMinimizer::Unchanged;
            if (w.state == Done || w.state == Fixed) continue;
            iret = w.improve(smallsteps,ytol);
            if (iret == OneDimMinimizer::Unchanged) {
                printf("\tMinimized %s:  Unchanged. NLL = %.8f\n", w.var().GetName(), nll_->getVal());
                w.state = Done;
                doneWorkers.push_front(&w);
            } else {
                printf("\tMinimized %s:  Changed. NLL = %.8f\n", w.var().GetName(), nll_->getVal());
                w.state = Active;
                newstate = Active;
            }
        }
        if (newstate == Done) {
            std::list<Worker*>::iterator it = doneWorkers.begin();
            while( it != doneWorkers.end()) {
                Worker &w = **it;
                OneDimMinimizer::ImproveRet iret = w.improve(smallsteps,ytol,0,/*force=*/true);
                if (iret == OneDimMinimizer::Unchanged) {
                    printf("\tMinimized %s:  Unchanged. NLL = %.8f\n", w.var().GetName(), nll_->getVal());
                    ++it;
                } else {
                    printf("\tMinimized %s:  Changed. NLL = %.8f\n", w.var().GetName(), nll_->getVal());
                    w.state = Active;
                    newstate = Active;
                    it = doneWorkers.erase(it);
                    break;
                }
            }
        }
        printf("End of loop. New state is %s\n",(newstate == Done ? "DONE" : "ACTIVE"));
        if (state_ == Done && newstate == Done) {
            //std::cout << "Converged after " << i << " big steps" << std::endl;
            return true;
        }
        state_ = newstate;
    }
    //std::cout << "Did not converge after " << bigsteps << " big steps" << std::endl;
    return false;
}


