FROM python:3

LABEL maintainer="kmaddela@nyit.edu" institution="New York Institute of Technology"

WORKDIR /usr/src/app
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
VOLUME /usr/src/app/output

CMD [ "bash" ]
