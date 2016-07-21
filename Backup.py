# for j in range(0, len(startTime)):
#     maxSS = 0
#     timeStamp = startTime[j]
#     # Get the index of the start time stamp
#     for k in range(previousK, len(csvContent)):
#         if csvContent[k][0] == timeStamp:
#             previousK = k
#             break
#     for i in range(self.timeFrame):
#         if k + i >= len(csvContent):
#             break
#         previousK = 0
#         sizeSum = 0
#         countSum = 0
#         ccnList.append(csvContent[k + i][-3])
#         if csvContent[k + i][1] > maxSS:
#             maxSS = csvContent[k + i][1]
#         for m in range(0, 20):
#             sizeSum += sizeList[m] * float(csvContent[k + i][binPos + m])
#             size += 0.5
#             countSum += float(csvContent[k + i][binPos + m])
#             # get the average size for each scan
#         if countSum == 0:
#             aveSizeList.append(0)
#         else:
#             aveSizeList.append(sizeSum / countSum)
#     self.ssList.append(maxSS)

# Get the rest of the csv data into the data so that there is extra data for later alignment
# extraCCNList = []
# extraAveSizeList = []
# for k in range(previousK, len(csvContent)):
#     if csvContent[k][0] == self.endTimeEntries[-1]:
#         previousK = k
#         break
# k += 1
# for i in range(self.timeFrame):
#     if k + i >= len(csvContent):
#         break
#     previousK = 0
#     sizeSum = 0
#     countSum = 0
#     extraCCNList.append(csvContent[k + i][-3])
#     for m in range(0, 20):
#         sizeSum += sizeList[m] * float(csvContent[k + i][binPos + m])
#         size += 0.5
#         countSum += float(csvContent[k + i][binPos + m])
#         # get the average size for each scan
#     if countSum == 0:
#         extraAveSizeList.append(0)
#     else:
#         extraAveSizeList.append(sizeSum / countSum)

# timeStamp = startTime[0]
# timeStamp = datetime.strptime(timeStamp, "%I:%M:%S")
# smpsCcnList = []
# for i in range(0, width):
#     for j in range(self.timeFrame):
#         aLine = [smpsList[j][1]] + [timeStamp.time()] + [float(smpsList[j][0])] + [
#             float(smpsList[j][i + 2])] + [float(ccnList[self.timeFrame * i + j])] + [float(aveSizeList[self.timeFrame * i + j])]
#         timeStamp = timeStamp + timedelta(seconds=1)
#         smpsCcnList.append(aLine)