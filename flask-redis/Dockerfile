FROM python:3.4-alpine
COPY requirements.txt .
RUN pip install -r requirements.txt \
&& rm requirements.txt
ADD code /code
WORKDIR /code
CMD ["python", "app.py"]