class Interval_Data
{
  public:
    double x;         /* original input */
    double lb;        /* lower bound */
    double ub;        /* upper bound */

    Interval_Data() {}
    ~Interval_Data() {}

    Interval_Data& operator=(const Interval_Data& rhs)
    {
        if(this == &rhs) return *this;

        this->x = rhs.x;
        this->lb = rhs.lb;
        this->ub = rhs.ub;

        return *this;
    }
};
