// Classification losses: L1 and Median. Extracted from cpf.cpp.
#include "cpf.hpp"

void ClassificationRPF::L1_loss(Split &split)
{
  split.min_sum = 0;
  split.M_s = split.sum_s / split.I_s.size();
  split.M_b = split.sum_b / split.I_b.size();
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_s[p]) - std::fabs((*split.Y)[individual][p]);
    for (auto individual : split.I_b)
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_b[p]) - std::fabs((*split.Y)[individual][p]);
  }
}

void ClassificationRPF::median_loss(Split &split)
{
  split.min_sum = 0;
  split.M_s = calcMedian(*split.Y, split.I_s);
  split.M_b = calcMedian(*split.Y, split.I_b);
  for (size_t p = 0; p < value_size; ++p)
  {
    for (auto individual : split.I_s)
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_s[p]) - std::fabs((*split.Y)[individual][p]);
    for (auto individual : split.I_b)
      split.min_sum += std::fabs((*split.Y)[individual][p] - split.M_b[p]) - std::fabs((*split.Y)[individual][p]);
  }
}


